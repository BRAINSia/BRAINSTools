#include "BRAINSMultiSTAPLECLP.h"
#include "itkIO.h"
#include "itkTransformFileReader.h"
#include "itkTransformFileWriter.h"
#include "itkIdentityTransform.h"
#include "itkResampleImageFilter.h"
#include "itkLabelImageGaussianInterpolateImageFunction.h"
#include "itkMultiLabelSTAPLEImageFilter.h"
#include "vnl/vnl_matlab_write.h"
#include <sstream>
#include <vector>

template <typename TImage>
void
printImageStats(const TImage *image)
{
  typename TImage::RegionType region = image->GetLargestPossibleRegion();
  std::cout << "size["
            << region.GetSize()[0] << ","
            << region.GetSize()[1] << ","
            << region.GetSize()[2] << "] Spacing["
            << image->GetSpacing()[0] << ","
            << image->GetSpacing()[1] << ","
            << image->GetSpacing()[2] << "] Origin ["
            << image->GetOrigin()[0] << ","
            << image->GetOrigin()[1] << ","
            << image->GetOrigin()[2] << "]" << std::endl;
}

int main(int argc, char *argv[])
{
  PARSE_ARGS;

  if( inputCompositeT1Volume == "" )
    {
    std::cerr << "Missing required Composite T1 Volume "
              << "use --inputCompositeT1Volume flag to specify"
              << std::endl;
    return 1;
    }
  if( inputLabelVolume.size() <= 0 )
    {
    std::cerr << "Missing input label volumes "
              << "use --inputLabelVolume <name> to add a label volume"
              << std::endl;
    return 1;
    }

  if( inputTransform.size() > 0 &&
      (inputLabelVolume.size() != inputTransform.size() ) )
    {
    std::cerr << "Transform list should have same number of" << std::endl
              << "members as as the input label volumes list" << std::endl;
    return 1;
    }

  if( outputMultiSTAPLE == "" )
    {
    std::cerr << "Missing outputMultiSTAPLE image file name" << std::endl;
    return 1;
    }

  typedef itk::Image<unsigned short, 3> USImageType;

  typedef std::vector<USImageType::Pointer> ImageList;
  ImageList inputLabelVolumes;
  for( std::vector<std::string>::const_iterator it = inputLabelVolume.begin();
       it != inputLabelVolume.end(); ++it )
    {
    USImageType::Pointer labelVolume;
    std::cout << "Reading " << (*it) << std::endl;
    try
      {
      labelVolume = itkUtil::ReadImage<USImageType>( (*it) );
      }
    catch( itk::ExceptionObject & err )
      {
      std::cerr << err << std::endl;
      return 1;
      }
    inputLabelVolumes.push_back(labelVolume);
    }

  ImageList transformedLabelVolumes;

  // resample all input label images into a common space defined by
  // the input Composite volume.
  if( skipResampling )
    {
    for( ImageList::const_iterator it = inputLabelVolumes.begin();
         it != inputLabelVolumes.end(); ++it )
      {
      transformedLabelVolumes.push_back( (*it) );
      }
    }
  else
    {
    USImageType::Pointer compositeVolume;
    try
      {
      std::cout << "Reading Composite Volume " << inputCompositeT1Volume
                << std::endl;
      compositeVolume = itkUtil::ReadImage<USImageType>(inputCompositeT1Volume);
      }
    catch( itk::ExceptionObject & err )
      {
      std::cerr << err << std::endl;
      return 1;
      }
    printImageStats<USImageType>(compositeVolume);

    typedef std::vector<itk::TransformFileReader::TransformPointer>
      TransformListType;

    TransformListType inputTransforms;

    if( inputTransform.size() > 0 )
      {
      for( std::vector<std::string>::const_iterator it = inputTransform.begin();
           it != inputTransform.end(); ++it )
        {
        itk::TransformFileReader::TransformPointer curTransform;

        itk::TransformFileReader::Pointer reader =
          itk::TransformFileReader::New();
        std::cout << "Reading " << (*it) << std::endl;
        reader->SetFileName( (*it) );
        try
          {
          reader->Update();
          }
        catch( itk::ExceptionObject & err )
          {
          std::cerr << err << std::endl;
          return 1;
          }
        curTransform = reader->GetTransformList()->front();
        inputTransforms.push_back(curTransform);
        }
      }
    else
      {
      std::cout << "No transforms specified, using Identity" << std::endl;
      // fake it with identity transforms
      typedef itk::IdentityTransform<double, 3> IDTransformType;

      IDTransformType::Pointer                   idXfrm = IDTransformType::New();
      itk::TransformFileReader::TransformPointer baseXfrm = idXfrm.GetPointer();
      for( std::vector<std::string>::const_iterator it = inputLabelVolume.begin();
           it != inputLabelVolume.end(); ++it )
        {
        inputTransforms.push_back(baseXfrm);
        }
      }
    // set up interpolator function
    // NOTE see ANTS/Examples/make_interpolator_snip.tmp line 113 --
    // the sigma defaults to the image spacing apparently, but the
    // sigma can also be specified on the command line.
    typedef std::less<itk::NumericTraits<unsigned char>::RealType> ucharLess;
    typedef itk::LabelImageGaussianInterpolateImageFunction<USImageType, double, ucharLess>
      InterpolationFunctionType;
    InterpolationFunctionType::Pointer interpolateFunc =
      InterpolationFunctionType::New();
    double                   sigma[3];
    USImageType::SpacingType spacing = compositeVolume->GetSpacing();
    for( unsigned i = 0; i < 3; ++i )
      {
      sigma[i] = spacing[i];
      }
    interpolateFunc->SetParameters(sigma, 4.0);

    std::vector<std::string>::const_iterator nameIt = inputLabelVolume.begin();
    TransformListType::const_iterator        xfrmIt = inputTransforms.begin();
    for( ImageList::const_iterator it = inputLabelVolumes.begin();
         it != inputLabelVolumes.end(); ++it, ++xfrmIt, ++nameIt )
      {
      USImageType::Pointer current = (*it);

      itk::TransformFileReader::TransformPointer curTransformBase = (*xfrmIt);

      typedef itk::ResampleImageFilter<USImageType, USImageType, double> ResampleFilterType;

      std::cout << "Resampling " << (*nameIt) << std::flush;

      const ResampleFilterType::TransformType *curTransform =
        dynamic_cast<const ResampleFilterType::TransformType *>(curTransformBase.GetPointer() );
      if( curTransform == 0 )
        {
        std::cerr << "Invalid transform " << curTransformBase << std::endl;
        exit(1);
        }
      ResampleFilterType::Pointer resampler = ResampleFilterType::New();
      try
        {
        resampler->SetInput(current);
        resampler->SetUseReferenceImage(true);
        resampler->SetReferenceImage(compositeVolume);
        resampler->SetInterpolator(interpolateFunc);
        resampler->SetTransform(curTransform);
        resampler->Update();
        }
      catch( itk::ExceptionObject & err )
        {
        std::cerr << err << std::endl;
        return 1;
        }
      std::cout << " done." << std::endl;
      if( resampledVolumePrefix != "" )
        {
        std::string namePart(itksys::SystemTools::GetFilenameName( (*nameIt) ) );
        std::string resampledName = resampledVolumePrefix;
        resampledName += namePart;
        std::cerr << "Writing " << resampledName << std::flush;
        try
          {
          itkUtil::WriteImage<USImageType>(resampler->GetOutput(), resampledName);
          }
        catch( itk::ExceptionObject & err )
          {
          std::cerr << err << std::endl;
          return 1;
          }
        std::cerr << " ... done." << std::endl;
        }
      printImageStats<USImageType>(resampler->GetOutput() );
      transformedLabelVolumes.push_back(resampler->GetOutput() );
      }
    }

  typedef itk::MultiLabelSTAPLEImageFilter<USImageType, USImageType> STAPLEFilterType;
  STAPLEFilterType::Pointer STAPLEFilter = STAPLEFilterType::New();
  STAPLEFilter->SetNumberOfThreads(1);

  if( labelForUndecidedPixels != -1 )
    {
    STAPLEFilter->SetLabelForUndecidedPixels(labelForUndecidedPixels);
    }
  for( ImageList::const_iterator it = transformedLabelVolumes.begin();
       it != transformedLabelVolumes.end(); ++it )
    {
    STAPLEFilter->PushBackInput( (*it) );
    }

  std::cout << "Running MultiLabel Staple filter " << std::flush;
  try
    {
    STAPLEFilter->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    std::cerr << err << std::endl;
    return 1;
    }
  USImageType::Pointer output = STAPLEFilter->GetOutput();

  std::cout << " done." << std::endl;

  try
    {
    std::cout << "Writing " << outputMultiSTAPLE << std::endl;
    itkUtil::WriteImage<USImageType>(output, outputMultiSTAPLE);
    }
  catch( itk::ExceptionObject & err )
    {
    std::cerr << err << std::endl;
    return 1;
    }

  if( outputConfusionMatrix != "" )
    {
    std::cout << "Writing " << outputConfusionMatrix << std::endl;
    std::ofstream out(outputConfusionMatrix.c_str() );
    if( !out.good() )
      {
      std::cerr << "Can't write Matlab confusion matrix file "
                << outputConfusionMatrix << std::endl;
      return 1;
      }
    for( unsigned int i = 0; i < transformedLabelVolumes.size(); ++i )
      {
      std::stringstream name;
      name << "confusionMat" << i;
      STAPLEFilterType::ConfusionMatrixType confusionMat =
        STAPLEFilter->GetConfusionMatrix(i);
      vnl_matlab_write(out, confusionMat.data_array(),
                       confusionMat.rows(),
                       confusionMat.cols(),
                       name.str().c_str() );
      }
    out.close();
    }
  return 0;
}
