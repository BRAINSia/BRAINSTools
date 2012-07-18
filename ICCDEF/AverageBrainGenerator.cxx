/*==================================================================

TODO:  NEED TO COMMENT WHAT THIS PROGRAM IS TO BE USED FOR

2009.08
Zhao,Yongqiang

==================================================================*/

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itksys/SystemTools.hxx"
#include "AverageBrainGeneratorCLP.h"
#include "itksys/Directory.hxx"
#include "itkAddImageFilter.h"
#include "itkMultiplyByConstantImageFilter.h"
#include "itkWarpImageFilter.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionIterator.h"
#include "itkIO.h"
#include "itkICCIterativeInverseDeformationFieldImageFilter.h"
#include "itkVectorIndexSelectionCastImageFilter.h"

int AverageBrainGenerator(int argc, char *argv[])
{
  PARSE_ARGS;

  const bool debug = true;

  if( debug )
    {
    std::cout << "=====================================================" << std::endl;
    std::cout << "Input Directory:     " <<  inputDirectory << std::endl;
    std::cout << "Template:            " <<  templateVolume << std::endl;
    std::cout << "Resolusion:          " <<  resolusion << std::endl;
    std::cout << "Iteration:           " <<  iteration << std::endl;
    std::cout << "Output Pixel Type:   " <<  pixelType << std::endl;
    std::cout << "Output Volume:       " <<  outputVolume << std::endl;
    std::cout << "=====================================================" << std::endl;
    }

  if( inputDirectory.size() == 0 )
    {
    std::cout << "Input Directory is misses!" << std::endl;
    exit(-1);
    }

  if( templateVolume.size() == 0 )
    {
    std::cout << "Template is missed!" << std::endl;
    exit(-1);
    }
  if( resolusion.size() == 0 )
    {
    std::cout << "Resolusion is missed!" << std::endl;
    exit(-1);
    }
  if( iteration.size() == 0 )
    {
    std::cout << "Iteration is missed!" << std::endl;
    exit(-1);
    }
  if( outputVolume.size() == 0 )
    {
    std::cout << "Output Volume is missed!" << std::endl;
    exit(-1);
    }

  const unsigned int Dimension = 3;
  typedef float                                   PixelType;
  typedef itk::Image<PixelType, Dimension>        ImageType;
  typedef itk::Vector<PixelType, Dimension>       VectorPixelType;
  typedef itk::Image<VectorPixelType,  Dimension> DeformationFieldType;
  DeformationFieldType::Pointer DeformationField = DeformationFieldType::New();

  ImageType::Pointer templateImage;
  templateImage = itkUtil::ReadImage<ImageType>(templateVolume);
  templateImage = itkUtil::OrientImage<ImageType>(templateImage,
                                                  itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_RAI);
  // Read Directory

  unsigned int numberOfFields = 0;

  DeformationField->CopyInformation(templateImage);
  DeformationField->SetRegions(templateImage->GetLargestPossibleRegion() );
  DeformationField->Allocate();

  DeformationFieldType::PixelType zeros;
  for( unsigned int j = 0; j < Dimension; j++ )
    {
    zeros[j] = 0.0;
    }

  itk::ImageRegionIterator<DeformationFieldType> it(DeformationField, DeformationField->GetRequestedRegion() );

  while( !it.IsAtEnd() )
    {
    it.Value() =  zeros;
    ++it;
    }

  std::cout << "Start..." << std::endl;
  itksys::Directory * dir = new itksys::Directory;
  std::string         subName = "iteration" + iteration + "resolution" + resolusion + "_forward.nii.gz";
  if( dir->Load(inputDirectory.c_str() ) )
    {
    std::cout << "Load..." << std::endl;
    for( unsigned long i = 0; i < dir->GetNumberOfFiles(); ++i )
      {
      std::string path = inputDirectory + "/" + dir->GetFile(i) + "/" + "forward";
      if( itksys::SystemTools::FileIsDirectory(path.c_str() ) )
        {
        // std::cout<<path<<std::endl;
        itksys::Directory * subDir = new itksys::Directory;
        if( subDir->Load(path.c_str() ) )
          {
          for( unsigned long j = 0; j < subDir->GetNumberOfFiles(); ++j )
            {
            if( itksys::SystemTools::StringEndsWith(subDir->GetFile(j), subName.c_str() ) )
              {
              // Compute the average displacement
              std::cout << subDir->GetFile(j) << std::endl;
              typedef itk::ImageFileReader<DeformationFieldType> DFReaderType;
              DFReaderType::Pointer df_Reader = DFReaderType::New();
              std::string           fileName = path + "/" + subDir->GetFile(j);
              df_Reader->SetFileName(fileName);
              df_Reader->Update();

              typedef itk::AddImageFilter<DeformationFieldType, DeformationFieldType,
                                          DeformationFieldType> AddImageType;
              AddImageType::Pointer adder = AddImageType::New();
              adder->SetInput1(DeformationField);
              adder->SetInput2(df_Reader->GetOutput() );
              adder->Update();
              DeformationField = adder->GetOutput();
              numberOfFields++;
              }
            }
          }
        }
      }
    }
  else
    {
    std::cout << "Can not open the directory!!!!" << std::endl;
    exit(-1);
    }

  if( numberOfFields < 3 )
    {
    std::cout << "NEED at least 3 data sets to make an average!" << std::endl;
    }

  typedef itk::MultiplyByConstantImageFilter<DeformationFieldType, float, DeformationFieldType> MultiplyImageType;
  MultiplyImageType::Pointer multi = MultiplyImageType::New();
  multi->SetInput(DeformationField);
  multi->SetConstant(1.0 / static_cast<float>(numberOfFields + 1) );
  multi->Update();

  // Compute the inverse of the average deformation field
  typedef itk::ICCIterativeInverseDeformationFieldImageFilter<DeformationFieldType,
                                                              DeformationFieldType> InverseDeformationFieldImageType;
  InverseDeformationFieldImageType::Pointer inverse = InverseDeformationFieldImageType::New();
  inverse->SetInput(multi->GetOutput() );
  inverse->SetStopValue(1.0e-6);
  inverse->SetNumberOfIterations(100);
  inverse->Update();

  // Write the displacement in each direction

  typedef itk::VectorIndexSelectionCastImageFilter<DeformationFieldType, ImageType> ComponentFilterType;
  ComponentFilterType::Pointer adaptor = ComponentFilterType::New();
  adaptor->SetInput(inverse->GetOutput() );

  std::string baseName = itksys::SystemTools::GetFilenameWithoutExtension(outputVolume.c_str() );
  char        ext[3][14] = { "_xdisp.nii.gz", "_ydisp.nii.gz", "_zdisp.nii.gz"};
  for( unsigned int extiter = 0; extiter < 3; extiter++ )
    {
    std::string componentFilename = baseName + ext[extiter];
    adaptor->SetIndex( extiter );
    adaptor->Update();

    typedef itk::ImageFileWriter<ImageType> ImageWriteType;
    ImageWriteType::Pointer writer = ImageWriteType::New();
    writer->SetInput(adaptor->GetOutput() );
    writer->SetFileName(componentFilename);
    writer->Update();
    }

  // Warp the templateImage with the avergae displacement
  typedef itk::WarpImageFilter<ImageType, ImageType, DeformationFieldType> WarpImageType;
  WarpImageType::Pointer warper = WarpImageType::New();

  warper->SetInput(templateImage);
  warper->SetOutputSpacing( DeformationField->GetSpacing() );
  warper->SetOutputOrigin( DeformationField->GetOrigin() );
  warper->SetOutputDirection( DeformationField->GetDirection() );
  warper->SetDeformationField( inverse->GetOutput() );

  warper->Update();

  if( pixelType == "uchar" )
    {
    typedef unsigned char                                 NewPixelType;
    typedef itk::Image<NewPixelType, Dimension>           NewImageType;
    typedef itk::CastImageFilter<ImageType, NewImageType> CastImageFilter;
    CastImageFilter::Pointer castFilter = CastImageFilter::New();
    castFilter->SetInput( warper->GetOutput() );
    castFilter->Update();

    typedef itk::ImageFileWriter<NewImageType> ImageWriteType;
    ImageWriteType::Pointer writer = ImageWriteType::New();
    writer->SetInput(castFilter->GetOutput() );
    writer->SetFileName(outputVolume);
    writer->Update();
    }

  else if( pixelType == "short" )
    {
    typedef short                                         NewPixelType;
    typedef itk::Image<NewPixelType, Dimension>           NewImageType;
    typedef itk::CastImageFilter<ImageType, NewImageType> CastImageFilter;
    CastImageFilter::Pointer castFilter = CastImageFilter::New();
    castFilter->SetInput( warper->GetOutput() );
    castFilter->Update();

    typedef itk::ImageFileWriter<NewImageType> ImageWriteType;
    ImageWriteType::Pointer writer = ImageWriteType::New();
    writer->SetInput(castFilter->GetOutput() );
    writer->SetFileName(outputVolume);
    writer->Update();
    }

  else if( pixelType == "ushort" )
    {
    typedef unsigned short                                NewPixelType;
    typedef itk::Image<NewPixelType, Dimension>           NewImageType;
    typedef itk::CastImageFilter<ImageType, NewImageType> CastImageFilter;
    CastImageFilter::Pointer castFilter = CastImageFilter::New();
    castFilter->SetInput( warper->GetOutput() );
    castFilter->Update();

    typedef itk::ImageFileWriter<NewImageType> ImageWriteType;
    ImageWriteType::Pointer writer = ImageWriteType::New();
    writer->SetInput(castFilter->GetOutput() );
    writer->SetFileName(outputVolume);
    writer->Update();
    }

  else if( pixelType == "int" )
    {
    typedef int                                           NewPixelType;
    typedef itk::Image<NewPixelType, Dimension>           NewImageType;
    typedef itk::CastImageFilter<ImageType, NewImageType> CastImageFilter;
    CastImageFilter::Pointer castFilter = CastImageFilter::New();
    castFilter->SetInput( warper->GetOutput() );
    castFilter->Update();

    typedef itk::ImageFileWriter<NewImageType> ImageWriteType;
    ImageWriteType::Pointer writer = ImageWriteType::New();
    writer->SetInput(castFilter->GetOutput() );
    writer->SetFileName(outputVolume);
    writer->Update();
    }

  else if( pixelType == "uint" )
    {
    typedef unsigned int                                  NewPixelType;
    typedef itk::Image<NewPixelType, Dimension>           NewImageType;
    typedef itk::CastImageFilter<ImageType, NewImageType> CastImageFilter;
    CastImageFilter::Pointer castFilter = CastImageFilter::New();
    castFilter->SetInput( warper->GetOutput() );
    castFilter->Update();

    typedef itk::ImageFileWriter<NewImageType> ImageWriteType;
    ImageWriteType::Pointer writer = ImageWriteType::New();
    writer->SetInput(castFilter->GetOutput() );
    writer->SetFileName(outputVolume);
    writer->Update();
    }

  else if( pixelType == "float" )
    {
    typedef itk::ImageFileWriter<ImageType> ImageWriteType;
    ImageWriteType::Pointer writer = ImageWriteType::New();
    writer->SetInput(warper->GetOutput() );
    writer->SetFileName(outputVolume);
    writer->Update();
    }

  else
    {
    std::cout << "ERROR:  Invalid pixelType" << std::endl;
    exit(-1);
    }

  return EXIT_SUCCESS;
}

int main( int argc, char *argv[] )
{
  // HACK:  BRAINS2 Masks are currently broken
  // The direction cosines are and the direction labels are not consistently being set.
  // itk::Brains2MaskImageIOFactory::RegisterOneFactory();

  return AverageBrainGenerator(argc, argv);
}
