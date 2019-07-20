/*=========================================================================
 *
 *  Copyright SINAPSE: Scalable Informatics for Neuroscience, Processing and Software Engineering
 *            The University of Iowa
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         http://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *=========================================================================*/
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
#include "itkMultiplyImageFilter.h"
#include "itkWarpImageFilter.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionIterator.h"
#include "itkIO.h"
#include "itkICCIterativeInverseDisplacementFieldImageFilter.h"
#include "itkVectorIndexSelectionCastImageFilter.h"
#include "BRAINSCommonLib.h"

int
AverageBrainGenerator( int argc, char * argv[] )
{
  PARSE_ARGS;
  BRAINSRegisterAlternateIO();

  const bool debug = true;

  if ( debug )
  {
    std::cout << "=====================================================" << std::endl;
    std::cout << "Input Directory:     " << inputDirectory << std::endl;
    std::cout << "Template:            " << templateVolume << std::endl;
    std::cout << "Resolusion:          " << resolusion << std::endl;
    std::cout << "Iteration:           " << iteration << std::endl;
    std::cout << "Output Pixel Type:   " << pixelType << std::endl;
    std::cout << "Output Volume:       " << outputVolume << std::endl;
    std::cout << "=====================================================" << std::endl;
  }

  if ( inputDirectory.size() == 0 )
  {
    std::cout << "Input Directory is misses!" << std::endl;
    exit( -1 );
  }

  if ( templateVolume.size() == 0 )
  {
    std::cout << "Template is missed!" << std::endl;
    exit( -1 );
  }
  if ( resolusion.size() == 0 )
  {
    std::cout << "Resolusion is missed!" << std::endl;
    exit( -1 );
  }
  if ( iteration.size() == 0 )
  {
    std::cout << "Iteration is missed!" << std::endl;
    exit( -1 );
  }
  if ( outputVolume.size() == 0 )
  {
    std::cout << "Output Volume is missed!" << std::endl;
    exit( -1 );
  }

  constexpr unsigned int Dimension = 3;
  using PixelType = float;
  using ImageType = itk::Image< PixelType, Dimension >;
  using VectorPixelType = itk::Vector< PixelType, Dimension >;
  using DisplacementFieldType = itk::Image< VectorPixelType, Dimension >;
  DisplacementFieldType::Pointer DisplacementField = DisplacementFieldType::New();

  ImageType::Pointer templateImage;
  templateImage = itkUtil::ReadImage< ImageType >( templateVolume );
  templateImage =
    itkUtil::OrientImage< ImageType >( templateImage, itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_RAI );
  // Read Directory

  unsigned int numberOfFields = 0;

  DisplacementField->CopyInformation( templateImage );
  DisplacementField->SetRegions( templateImage->GetLargestPossibleRegion() );
  DisplacementField->Allocate();

  DisplacementFieldType::PixelType zeros;
  for ( unsigned int j = 0; j < Dimension; j++ )
  {
    zeros[j] = 0.0;
  }

  itk::ImageRegionIterator< DisplacementFieldType > it( DisplacementField, DisplacementField->GetRequestedRegion() );

  while ( !it.IsAtEnd() )
  {
    it.Value() = zeros;
    ++it;
  }

  std::cout << "Start..." << std::endl;
  itksys::Directory * dir = new itksys::Directory;
  std::string         subName = "iteration" + iteration + "resolution" + resolusion + "_forward.nii.gz";
  if ( dir->Load( inputDirectory.c_str() ) )
  {
    std::cout << "Load..." << std::endl;
    for ( unsigned long i = 0; i < dir->GetNumberOfFiles(); ++i )
    {
      std::string path = inputDirectory + "/" + dir->GetFile( i ) + "/" + "forward";
      if ( itksys::SystemTools::FileIsDirectory( path.c_str() ) )
      {
        // std::cout<<path<<std::endl;
        itksys::Directory * subDir = new itksys::Directory;
        if ( subDir->Load( path.c_str() ) )
        {
          for ( unsigned long j = 0; j < subDir->GetNumberOfFiles(); ++j )
          {
            if ( itksys::SystemTools::StringEndsWith( subDir->GetFile( j ), subName.c_str() ) )
            {
              // Compute the average displacement
              std::cout << subDir->GetFile( j ) << std::endl;
              using DFReaderType = itk::ImageFileReader< DisplacementFieldType >;
              DFReaderType::Pointer df_Reader = DFReaderType::New();
              std::string           fileName = path + "/" + subDir->GetFile( j );
              df_Reader->SetFileName( fileName );
              df_Reader->Update();

              using AddImageType =
                itk::AddImageFilter< DisplacementFieldType, DisplacementFieldType, DisplacementFieldType >;
              AddImageType::Pointer adder = AddImageType::New();
              adder->SetInput1( DisplacementField );
              adder->SetInput2( df_Reader->GetOutput() );
              adder->Update();
              DisplacementField = adder->GetOutput();
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
    exit( -1 );
  }

  if ( numberOfFields < 3 )
  {
    std::cout << "NEED at least 3 data sets to make an average!" << std::endl;
  }

  using MultiplyImageType =
    itk::MultiplyImageFilter< DisplacementFieldType, itk::Image< float, Dimension >, DisplacementFieldType >;
  MultiplyImageType::Pointer multi = MultiplyImageType::New();
  multi->SetInput( DisplacementField );
  multi->SetConstant( 1.0 / static_cast< float >( numberOfFields + 1 ) );
  multi->Update();

  // Compute the inverse of the average deformation field
  using InverseDisplacementFieldImageType =
    itk::ICCIterativeInverseDisplacementFieldImageFilter< DisplacementFieldType, DisplacementFieldType >;
  InverseDisplacementFieldImageType::Pointer inverse = InverseDisplacementFieldImageType::New();
  inverse->SetInput( multi->GetOutput() );
  inverse->SetStopValue( 1.0e-6 );
  inverse->SetNumberOfIterations( 100 );
  inverse->Update();

  // Write the displacement in each direction

  using ComponentFilterType = itk::VectorIndexSelectionCastImageFilter< DisplacementFieldType, ImageType >;
  ComponentFilterType::Pointer adaptor = ComponentFilterType::New();
  adaptor->SetInput( inverse->GetOutput() );

  std::string baseName = itksys::SystemTools::GetFilenameWithoutExtension( outputVolume.c_str() );
  char        ext[3][14] = { "_xdisp.nii.gz", "_ydisp.nii.gz", "_zdisp.nii.gz" };
  for ( unsigned int extiter = 0; extiter < 3; extiter++ )
  {
    std::string componentFilename = baseName + ext[extiter];
    adaptor->SetIndex( extiter );
    adaptor->Update();

    using ImageWriteType = itk::ImageFileWriter< ImageType >;
    ImageWriteType::Pointer writer = ImageWriteType::New();
    writer->SetInput( adaptor->GetOutput() );
    writer->SetFileName( componentFilename );
    writer->Update();
  }

  // Warp the templateImage with the avergae displacement
  using WarpImageType = itk::WarpImageFilter< ImageType, ImageType, DisplacementFieldType >;
  WarpImageType::Pointer warper = WarpImageType::New();

  warper->SetInput( templateImage );
  warper->SetOutputSpacing( DisplacementField->GetSpacing() );
  warper->SetOutputOrigin( DisplacementField->GetOrigin() );
  warper->SetOutputDirection( DisplacementField->GetDirection() );
  warper->SetDisplacementField( inverse->GetOutput() );

  warper->Update();

  if ( pixelType == "uchar" )
  {
    using NewPixelType = unsigned char;
    using NewImageType = itk::Image< NewPixelType, Dimension >;
    using CastImageFilter = itk::CastImageFilter< ImageType, NewImageType >;
    CastImageFilter::Pointer castFilter = CastImageFilter::New();
    castFilter->SetInput( warper->GetOutput() );
    castFilter->Update();

    using ImageWriteType = itk::ImageFileWriter< NewImageType >;
    ImageWriteType::Pointer writer = ImageWriteType::New();
    writer->SetInput( castFilter->GetOutput() );
    writer->SetFileName( outputVolume );
    writer->Update();
  }

  else if ( pixelType == "short" )
  {
    using NewPixelType = short;
    using NewImageType = itk::Image< NewPixelType, Dimension >;
    using CastImageFilter = itk::CastImageFilter< ImageType, NewImageType >;
    CastImageFilter::Pointer castFilter = CastImageFilter::New();
    castFilter->SetInput( warper->GetOutput() );
    castFilter->Update();

    using ImageWriteType = itk::ImageFileWriter< NewImageType >;
    ImageWriteType::Pointer writer = ImageWriteType::New();
    writer->SetInput( castFilter->GetOutput() );
    writer->SetFileName( outputVolume );
    writer->Update();
  }

  else if ( pixelType == "ushort" )
  {
    using NewPixelType = unsigned short;
    using NewImageType = itk::Image< NewPixelType, Dimension >;
    using CastImageFilter = itk::CastImageFilter< ImageType, NewImageType >;
    CastImageFilter::Pointer castFilter = CastImageFilter::New();
    castFilter->SetInput( warper->GetOutput() );
    castFilter->Update();

    using ImageWriteType = itk::ImageFileWriter< NewImageType >;
    ImageWriteType::Pointer writer = ImageWriteType::New();
    writer->SetInput( castFilter->GetOutput() );
    writer->SetFileName( outputVolume );
    writer->Update();
  }

  else if ( pixelType == "int" )
  {
    using NewPixelType = int;
    using NewImageType = itk::Image< NewPixelType, Dimension >;
    using CastImageFilter = itk::CastImageFilter< ImageType, NewImageType >;
    CastImageFilter::Pointer castFilter = CastImageFilter::New();
    castFilter->SetInput( warper->GetOutput() );
    castFilter->Update();

    using ImageWriteType = itk::ImageFileWriter< NewImageType >;
    ImageWriteType::Pointer writer = ImageWriteType::New();
    writer->SetInput( castFilter->GetOutput() );
    writer->SetFileName( outputVolume );
    writer->Update();
  }

  else if ( pixelType == "uint" )
  {
    using NewPixelType = unsigned int;
    using NewImageType = itk::Image< NewPixelType, Dimension >;
    using CastImageFilter = itk::CastImageFilter< ImageType, NewImageType >;
    CastImageFilter::Pointer castFilter = CastImageFilter::New();
    castFilter->SetInput( warper->GetOutput() );
    castFilter->Update();

    using ImageWriteType = itk::ImageFileWriter< NewImageType >;
    ImageWriteType::Pointer writer = ImageWriteType::New();
    writer->SetInput( castFilter->GetOutput() );
    writer->SetFileName( outputVolume );
    writer->Update();
  }

  else if ( pixelType == "float" )
  {
    using ImageWriteType = itk::ImageFileWriter< ImageType >;
    ImageWriteType::Pointer writer = ImageWriteType::New();
    writer->SetInput( warper->GetOutput() );
    writer->SetFileName( outputVolume );
    writer->Update();
  }

  else
  {
    std::cout << "ERROR:  Invalid pixelType" << std::endl;
    exit( -1 );
  }

  return EXIT_SUCCESS;
}

int
main( int argc, char * argv[] )
{
  // HACK:  BRAINS2 Masks are currently broken
  // The direction cosines are and the direction labels are not consistently being set.
  // itk::Brains2MaskImageIOFactory::RegisterOneFactory();

  return AverageBrainGenerator( argc, argv );
}
