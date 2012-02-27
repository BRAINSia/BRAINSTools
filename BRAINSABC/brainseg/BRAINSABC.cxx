#include "itkOutputWindow.h"
#include "itkTextOutput.h"
#include "itkTimeProbe.h"
#include <exception>
#include <iostream>
#include <string>
#include <fstream>
#include <map>
#include <vector>
#include <algorithm>

#include "itkNormalizedCorrelationImageToImageMetric.h"
#include "itkCastImageFilter.h"
#include "itkImage.h"
#include "itkImageFileWriter.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkWindowedSincInterpolateImageFunction.h"
#include "itkNumericTraits.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkVersion.h"
#include "itksys/SystemTools.hxx"
#include "PrettyPrintTable.h"
#include "DenoiseFiltering.h"
#include "itkImageDuplicator.h"

typedef itk::Image<float, 3>         FloatImageType;
typedef itk::Image<unsigned char, 3> ByteImageType;
typedef itk::Image<short, 3>         ShortImageType;

typedef FloatImageType::Pointer FloatImagePointer;
typedef ByteImageType::Pointer  ByteImagePointer;
typedef ShortImageType::Pointer ShortImagePointer;

#include "mu.h"
#include "EMSParameters.h"
#include "AtlasDefinition.h"
#include <vector>
#include "Log.h"

#include <itksys/SystemTools.hxx>

#include <iostream>
#include <string>
#include <sstream>
#include <map>

#include <stdlib.h>

#include <StandardizeMaskIntensity.h>
#include "BRAINSABCCLP.h"
#include "BRAINSThreadControl.h"

// Use manually instantiated classes for the big program chunks
#define MU_MANUAL_INSTANTIATION
#include "EMSegmentationFilter.h"
#include "AtlasRegistrationMethod.h"
#undef MU_MANUAL_INSTANTIATION

// From an inputVector create an outputVector that only contains the unique
// elements.
template <class T>
T RemoveDuplicates(T & origVect, const std::vector<bool> & dupsList)
{
  T outVector(0);

  outVector.reserve(origVect.size() );
  for( size_t inlist = 0; inlist < origVect.size(); inlist++ )
    {
    if( dupsList[inlist] == false )
      {
      outVector.push_back(origVect[inlist]);
      }
    }
  return outVector;
}

//
// Get a version of the filename that does not include the preceeding path, or
// the image file extensions.
static std::string GetStripedImageFileNameExtension(const std::string & ImageFileName)
{
  std::vector<std::string> ExtensionsToRemove(6);

  ExtensionsToRemove[0] = ".gz";
  ExtensionsToRemove[1] = ".nii";
  ExtensionsToRemove[2] = ".hdr";
  ExtensionsToRemove[3] = ".img";
  ExtensionsToRemove[4] = ".gipl";
  ExtensionsToRemove[5] = ".nrrd";

  std::string returnString = itksys::SystemTools::GetFilenameName( ImageFileName );
  for( size_t s = 0; s < ExtensionsToRemove.size(); s++ )
    {
    size_t rfind_location = returnString.rfind(ExtensionsToRemove[s]);
    if( ( rfind_location != std::string::npos )
        && ( rfind_location == ( returnString.size() - ExtensionsToRemove[s].size() ) )
        )
      {
      returnString.replace(rfind_location, ExtensionsToRemove[s].size(), "");
      }
    }
  return returnString;
}

static FloatImageType::Pointer CopyOutputImage(FloatImageType::Pointer img )
{
  // muLogMacro(<< "CopyOutputImage" << std::endl );

  FloatImageType::Pointer outimg = FloatImageType::New();

  outimg->CopyInformation(img);
  outimg->SetRegions( img->GetLargestPossibleRegion() );
  outimg->Allocate();

  typedef itk::ImageRegionIterator<FloatImageType> InternalIteratorType;
  InternalIteratorType inputIter( img, img->GetLargestPossibleRegion() );

  typedef itk::ImageRegionIterator<FloatImageType> OutputIteratorType;
  OutputIteratorType outputIter( outimg, outimg->GetLargestPossibleRegion() );

  inputIter.GoToBegin();
  outputIter.GoToBegin();

  while( !inputIter.IsAtEnd() )
    {
    outputIter.Set( static_cast<FloatImageType::PixelType>( inputIter.Get() ) );
    ++inputIter;
    ++outputIter;
    }

  return outimg;
}

// /////
//
//
//
static std::vector<FloatImageType::Pointer> ResampleImageList(
  const std::string & resamplerInterpolatorType,
  const std::vector<FloatImageType::Pointer> & inputImageList,
  const std::vector<GenericTransformType::Pointer> & intraSubjectTransforms)
{
  // Clear image list
  std::vector<FloatImageType::Pointer> outputImageList;

  outputImageList.clear();
  outputImageList.resize( inputImageList.size() );
  outputImageList[0] = CopyOutputImage(inputImageList[0]);

  const FloatImageType::PixelType outsideFOVCode = vnl_huge_val( static_cast<FloatImageType::PixelType>( 1.0f ) );
  // Resample the other images
  for( unsigned int i = 1; i < inputImageList.size(); i++ )
    {
    muLogMacro(<< "Resampling input image " << i + 1 << "." << std::endl);
    typedef itk::ResampleImageFilter<FloatImageType, FloatImageType> ResampleType;
    typedef ResampleType::Pointer                                    ResamplePointer;
    ResamplePointer resampler = ResampleType::New();
    resampler->SetInput(inputImageList[i]);
    resampler->SetTransform(intraSubjectTransforms[i]);

    if( resamplerInterpolatorType == "BSpline" )
      {
      typedef itk::BSplineInterpolateImageFunction<FloatImageType, double, double>
        SplineInterpolatorType;

      // Spline interpolation, only available for input images, not
      // atlas
      SplineInterpolatorType::Pointer splineInt
        = SplineInterpolatorType::New();
      splineInt->SetSplineOrder(5);
      resampler->SetInterpolator(splineInt);
      }
    else if( resamplerInterpolatorType == "WindowedSinc" )
      {
      typedef itk::ConstantBoundaryCondition<FloatImageType>
        BoundaryConditionType;
      static const unsigned int WindowedSincHammingWindowRadius = 5;
      typedef itk::Function::HammingWindowFunction<
          WindowedSincHammingWindowRadius, double, double> WindowFunctionType;
      typedef itk::WindowedSincInterpolateImageFunction
        <FloatImageType,
         WindowedSincHammingWindowRadius,
         WindowFunctionType,
         BoundaryConditionType,
         double>    WindowedSincInterpolatorType;
      WindowedSincInterpolatorType::Pointer windowInt
        = WindowedSincInterpolatorType::New();
      resampler->SetInterpolator(windowInt);
      }
    else // Default to m_UseNonLinearInterpolation == "Linear"
      {
      typedef itk::LinearInterpolateImageFunction<FloatImageType, double>
        LinearInterpolatorType;
      LinearInterpolatorType::Pointer linearInt
        = LinearInterpolatorType::New();
      resampler->SetInterpolator(linearInt);
      }

    resampler->SetDefaultPixelValue(outsideFOVCode);
    resampler->SetOutputParametersFromImage(inputImageList[0]);
    resampler->Update();

    // Zero the mask region outside FOV and also the intensities with
    // outside
    // FOV code
    typedef itk::ImageRegionIterator<FloatImageType> InternalIteratorType;

    FloatImageType::Pointer tmp = resampler->GetOutput();
    InternalIteratorType    tmpIt( tmp, tmp->GetLargestPossibleRegion() );

    // HACK:  We can probably remove the mask generation from here.
    // The FOV mask, regions where intensities in all channels do not
    // match FOV code
    ByteImageType::Pointer intraSubjectFOVIntersectionMask = NULL;
    intraSubjectFOVIntersectionMask = ByteImageType::New();
    intraSubjectFOVIntersectionMask->CopyInformation(inputImageList[0]);
    intraSubjectFOVIntersectionMask->SetRegions( inputImageList[0]->GetLargestPossibleRegion() );
    intraSubjectFOVIntersectionMask->Allocate();
    intraSubjectFOVIntersectionMask->FillBuffer(1);
    typedef itk::ImageRegionIterator<ByteImageType> MaskIteratorType;
    MaskIteratorType maskIt( intraSubjectFOVIntersectionMask,
                             intraSubjectFOVIntersectionMask->GetLargestPossibleRegion() );
    maskIt.GoToBegin();
    tmpIt.GoToBegin();
    while( !maskIt.IsAtEnd() )
      {
      if( tmpIt.Get() == outsideFOVCode )  // Voxel came from outside
      // the original FOV during
      // registration, so
      // invalidate it.
        {
        maskIt.Set(0); // Set it as an invalid voxel in
        // intraSubjectFOVIntersectionMask
        tmpIt.Set(0);  // Set image intensity value to zero.
        }
      ++maskIt;
      ++tmpIt;
      }

    // Add the image
    outputImageList[i] = CopyOutputImage(tmp);
    }
  return outputImageList;
}

static void RescaleFunctionLocal( std::vector<FloatImageType::Pointer> & localList)
{
  for( unsigned int i = 0; i < localList.size(); i++ )
    {
    typedef itk::RescaleIntensityImageFilter<FloatImageType, FloatImageType>
      RescaleType;
    RescaleType::Pointer rescaler = RescaleType::New();
    rescaler->SetOutputMinimum(1);
    rescaler->SetOutputMaximum(MAX_IMAGE_OUTPUT_VALUE);
#define INPLACE_RESCALER 0 // HACK Test this out
#if defined( INPLACE_RESCALER )
    rescaler->SetInPlace(true);
#endif
    FloatImageType::Pointer tmp = localList[i];
    rescaler->SetInput(tmp);
    rescaler->Update();
#if defined( INPLACE_RESCALER )
    localList[i] = rescaler->GetOutput();
#else
    FloatImageType::SizeType size = localList[0]->GetLargestPossibleRegion().GetSize();

    FloatImagePointer         rImg = rescaler->GetOutput();
    FloatImageType::IndexType ind;
    // T.O.D.O.  This could be done in-place using the -inplace flag of the
    // rescaleImageIntensityFilter.
      {
#pragma omp parallel for
      for( long kk = 0; kk < (long)size[2]; kk++ )
        {
        for( long jj = 0; jj < (long)size[1]; jj++ )
          {
          for( long ii = 0; ii < (long)size[0]; ii++ )
            {
            const ProbabilityImageIndexType currIndex = {{ii, jj, kk}};
            tmp->SetPixel( currIndex, rImg->GetPixel(ind) );
            }
          }
        }
      }
#endif
    }
}

static std::vector<bool> FindDuplicateImages(const std::vector<FloatImagePointer> candidateSameImageList )
{
  const double THREASHOLD_CUTOFF = 0.999;  // Images with higher correlation are

  // considered duplicate.

  typedef itk::IdentityTransform<double, FloatImageType::ImageDimension> IDTYPE;
  IDTYPE::Pointer myID = IDTYPE::New();
  typedef itk::LinearInterpolateImageFunction<FloatImageType, double> InterpType;
  InterpType::Pointer myInterp = InterpType::New();

  std::vector<bool> isDuplicated(candidateSameImageList.size(), false);
  for( unsigned int start = 0; start < candidateSameImageList.size(); start++ )
    {
    for( unsigned int q = start + 1; q < candidateSameImageList.size(); q++ )
      {
      typedef itk::NormalizedCorrelationImageToImageMetric<FloatImageType, FloatImageType> NormalizerType;
      NormalizerType::Pointer myNormalizer = NormalizerType::New();
      myNormalizer->SetFixedImage(candidateSameImageList[start]);
      myNormalizer->SetMovingImage(candidateSameImageList[q]);
      myNormalizer->SetTransform(myID);
      myNormalizer->SetFixedImageRegion( candidateSameImageList[start]->GetBufferedRegion() );
      myInterp->SetInputImage(candidateSameImageList[q]);
      myNormalizer->SetInterpolator(myInterp);
      myNormalizer->Initialize();
      const double correlationValue = vcl_abs(myNormalizer->GetValue(myID->GetParameters() ) );
      std::cout << "Correlation value between image " << start
                << " and image " << q << ": " << correlationValue << std::endl;
      if( correlationValue > THREASHOLD_CUTOFF )
        {
        isDuplicated[q] = true;
        }
      }
    }
  for( unsigned int q = 0; q < candidateSameImageList.size(); q++ )
    {
    std::cout << "Removing highly correlated image " << static_cast<int>(isDuplicated[q]) << std::endl;
    }
  return isDuplicated;
}

class EmptyVectorException
{
public:
  EmptyVectorException(const char* pStr = "The list of input images was empty.  Nothing to averge.") :
    pMessage(pStr)
  {
  }

  const char * what() const
  {
    return pMessage;
  }

private:
  const char * pMessage;
};

// Take a list of coregistered images, all of the same type (T1,T2) and return the average image.
static FloatImageType::Pointer AverageImageList(
  const std::vector<FloatImageType::Pointer> & inputImageList)
{
  if( inputImageList.size() == 0 )
    {
    // No images, something went wrong.
    throw EmptyVectorException();
    }
  if( inputImageList.size() == 1 )
    {
    // Only one image, nothing to average.
    return inputImageList[0];
    }

  // Create an image iterator over the first image.  Use that iterator to get an index into the other
  // images, sum each of the voxel values and divide by the number of input images and set the output
  // voxel at this index to that value.

  // Duplicate the first input image to use as an output image.
  typedef itk::ImageDuplicator<FloatImageType> DuplicatorType;
  DuplicatorType::Pointer duplicator = DuplicatorType::New();
  duplicator->SetInputImage(inputImageList[0]);
  duplicator->Update();
  FloatImageType::Pointer averageImage = duplicator->GetOutput();

  // Create an image iterator over the first image.
  typedef itk::ImageRegionIterator<FloatImageType> ImageRegionIteratorType;
  ImageRegionIteratorType imgItr( inputImageList[0], inputImageList[0]->GetRequestedRegion() );
  // Loop over the voxels calculating the averages.
  for( imgItr.GoToBegin(); !imgItr.IsAtEnd(); ++imgItr )
    {
    const FloatImageType::IndexType & idx = imgItr.GetIndex();
    FloatImageType::PixelType         avgValue = 0;
    const unsigned int                inputListSize = inputImageList.size();
    const float                       invListSize = 1.0F / static_cast<float>(inputListSize);
    for( unsigned int j = 0; j < inputImageList.size(); ++j )
      {
      avgValue += inputImageList[j]->GetPixel(idx);
      }
    averageImage->SetPixel(idx, static_cast<FloatImageType::PixelType>(avgValue * invListSize) );
    }

  return averageImage;
}

int main(int argc, char * *argv)
{
  PARSE_ARGS;
  const BRAINSUtils::StackPushITKDefaultNumberOfThreads TempDefaultNumberOfThreadsHolder(numberOfThreads);

  // TODO:  Need to figure out how to conserve memory better during the running
  // of this application:  itk::DataObject::GlobalReleaseDataFlagOn();
  itk::OutputWindow::SetInstance( itk::TextOutput::New() );

  typedef EMSegmentationFilter<FloatImageType, FloatImageType> SegFilterType;

  // Check the parameters for valid values
  bool AllSimpleParameterChecksValid = true;
  if( maxIterations < 1 )
    {
    muLogMacro( <<  "Warning: "
                << "--maxIterations set to 0, so only initialization with priors will be completed." << std::endl );
    }
  if( inputVolumes.size() == 0 )
    {
    muLogMacro( <<  "ERROR: "
                << "Must specify --inputVolumes" << std::endl );
    AllSimpleParameterChecksValid = false;
    }
  if( outputVolumes.size() != 1 && outputVolumes.size() != inputVolumes.size() )
    {
    std::cerr << inputVolumes.size() << " images in input volumeslist, but "
              << outputVolumes.size() << " names in output volumes list"
              << "OR it must be exactly 1, and be the template for writing files."
              << std::endl;
    return EXIT_FAILURE;
    }
  if( inputVolumeTypes.size() != inputVolumes.size() )
    {
    muLogMacro( <<  "ERROR: "
                << "--inputVolumeTypes and --inputVolumes must"
                << " have the same number of elements" << std::endl );
    AllSimpleParameterChecksValid = false;
    }
  if( atlasDefinition == "" )
    {
    muLogMacro( <<  "Error: "
                << "--atlasDefinition <xml atlas def> required"
                << std::endl );
    AllSimpleParameterChecksValid = false;
    }
  if( outputDir == "" )
    {
    muLogMacro( <<  "ERROR: "
                << "outputDir must be specified" << std::endl );
    AllSimpleParameterChecksValid = false;
    }
  if( AllSimpleParameterChecksValid == false )
    {
    muLogMacro( << "ERROR:  Commanline arguments are not valid." << std::endl );
    GENERATE_ECHOARGS;
    return EXIT_FAILURE;
    }
  AtlasDefinition atlasDefinitionParser;
  try
    {
    atlasDefinitionParser.InitFromXML(atlasDefinition);
    }
  catch( ... )
    {
    muLogMacro( <<  "Error reading Atlas Definition from "
                << atlasDefinition
                << std::endl );
    return EXIT_FAILURE;
    }
  ;
  atlasDefinitionParser.DebugPrint();
  // Create and start a new timer (for the whole process)
  //  EMSTimer* timer = new EMSTimer();
  itk::TimeProbe timer;
  timer.Start();

  // Directory separator string
  std::string separator = std::string("/");
  //  separator[0] = MU_DIR_SEPARATOR; // always a '/' -- windows can handle
  // that fine...

  // Make sure last character in output directory string is a separator
  if( outputDir[outputDir.size() - 1] != '/' /* MU_DIR_SEPARATOR */ )
    {
    outputDir += separator;
    }

  // Create the output directory, stop if it does not exist
  // if(!mu::create_dir(outputDir.c_str()))
  if( !itksys::SystemTools::MakeDirectory( outputDir.c_str() ) )
    {
    muLogMacro( << "ERROR: Could not create requested output directory " << outputDir << std::endl );
    return EXIT_FAILURE;
    }

  // Set up the logger
    {
    const std::string logfn = outputDir + defaultSuffix + ".log";
    ( mu::Log::GetInstance() )->EchoOn();
    ( mu::Log::GetInstance() )->SetOutputFileName( logfn.c_str() );
    }

  // Set up suffix string for images
  std::string fmt = outputFormat;
  std::string outext = ".mha";
  if( itksys::SystemTools::Strucmp(fmt.c_str(), "Nrrd") == 0 )
    {
    outext = ".nrrd";
    }
  else if( itksys::SystemTools::Strucmp(fmt.c_str(), "Meta") == 0 )
    {
    outext = ".mha";
    }
  else if( itksys::SystemTools::Strucmp(fmt.c_str(), "NIFTI") == 0 )
    {
    outext = ".nii.gz";
    }
  else
    {
    muLogMacro(<< "WARNING: output format unrecognized, using Meta format\n");
    }
  const std::string suffstr
    = std::string("_") + std::string(defaultSuffix) + outext;

  muLogMacro(<< "mu::brainseg\n");
  muLogMacro(<< "========================================\n");
  muLogMacro(<< "Program compiled on: " << __DATE__ << std::endl );
  muLogMacro(<< std::endl );

  muLogMacro(
    << "Hans J. Johnson - hans-johnson@uiowa.edu, has made significant"
    << " edits to this front end of the BRAINSABC system.\n");
  muLogMacro(
    << "Original application was written by Marcel Prastawa - "
    << "prastawa@sci.utah.edu, and is maintained as a separate program.\n");
  muLogMacro(<< "This software is provided for research purposes only\n");
  muLogMacro(<< std::endl );

  muLogMacro(<< "Using ITK version "
             << itk::Version::GetITKMajorVersion() << "."
             << itk::Version::GetITKMinorVersion() << "."
             << itk::Version::GetITKBuildVersion() << std::endl );
  muLogMacro(<< std::endl );

  // Write input parameters
  muLogMacro(<< "=== Parameters ===\n");
  muLogMacro(<< "Suffix: " << defaultSuffix << std::endl );
  muLogMacro(<< "Output Directory: " << outputDir << std::endl );
  muLogMacro(<< "Output Format: " << outputFormat << std::endl );
  muLogMacro(<< "Input images: \n");
  muLogMacro(
    << "Non-linear filtering, method: " << filterMethod << ", "
    << filterIteration
    << " iterations, dt = " << filterTimeStep << std::endl );

  const AtlasDefinition::TissueTypeVector & PriorNames
    = atlasDefinitionParser.TissueTypes();

  PrettyPrintTable AtlasDefTable;
  AtlasDefTable.add(0, 0, "Prior Names");
  AtlasDefTable.add(0, 1, ": [");
  AtlasDefTable.add(0, PriorNames.size() + 2 + 1, "]");
  SegFilterType::VectorType priorsWeightList;
  priorsWeightList.set_size( PriorNames.size() );
  AtlasDefTable.add(1, 0, "Prior weight scales");
  AtlasDefTable.add(1, 1, ": [");
  AtlasDefTable.add(1, PriorNames.size() + 2 + 1, "]");

  unsigned int currentRow = 0;
  for( unsigned int pwi = 0; pwi < PriorNames.size(); pwi++ )
    {
    AtlasDefTable.add(currentRow, 2 + pwi, PriorNames[pwi]);
    }

  currentRow++;
  for( unsigned int pwi = 0; pwi < PriorNames.size(); pwi++ )
    {
    priorsWeightList[pwi] = atlasDefinitionParser.GetWeight(PriorNames[pwi]);
    AtlasDefTable.add(currentRow, 2 + pwi, priorsWeightList[pwi], "%4.2f");
    }

  currentRow++;
  SegFilterType::IntVectorType priorLabelCodeVector;
  priorLabelCodeVector.set_size( PriorNames.size() );
  AtlasDefTable.add(currentRow, 0, "Prior Label Codes");
  AtlasDefTable.add(currentRow, 1, ": [");
  AtlasDefTable.add(currentRow, PriorNames.size() + 2 + 1, "]");
  for( unsigned int pwi = 0; pwi < PriorNames.size(); pwi++ )
    {
    priorLabelCodeVector[pwi] = atlasDefinitionParser.GetLabelCode(PriorNames[pwi]);
    AtlasDefTable.add(currentRow, 2 + pwi, priorLabelCodeVector[pwi], "%d");
    }

  currentRow++;
  SegFilterType::BoolVectorType priorIsForegroundPriorVector;
  priorIsForegroundPriorVector.resize( PriorNames.size() );
  AtlasDefTable.add(currentRow, 0, "Prior IsForeground");
  AtlasDefTable.add(currentRow, 1, ": [");
  AtlasDefTable.add(currentRow, PriorNames.size() + 2 + 1, "]");
  for( unsigned int pwi = 0; pwi < PriorNames.size(); pwi++ )
    {
    priorIsForegroundPriorVector[pwi] = atlasDefinitionParser.GetIsForegroundPrior(PriorNames[pwi]);
    AtlasDefTable.add(currentRow, 2 + pwi, priorIsForegroundPriorVector[pwi], "%d");
    }

  currentRow++;
  SegFilterType::IntVectorType priorGaussianClusterCountVector;
  priorGaussianClusterCountVector.set_size( PriorNames.size() );
  AtlasDefTable.add(currentRow, 0, "Prior Clusters");
  AtlasDefTable.add(currentRow, 1, ": [");
  AtlasDefTable.add(currentRow, PriorNames.size() + 2 + 1, "]");
  for( unsigned int pwi = 0; pwi < PriorNames.size(); pwi++ )
    {
    priorGaussianClusterCountVector[pwi] = atlasDefinitionParser.GetGaussianClusterCount(PriorNames[pwi]);
    AtlasDefTable.add(currentRow, 2 + pwi, priorGaussianClusterCountVector[pwi], "%d");
    }

  currentRow++;
  SegFilterType::BoolVectorType priorUseForBiasVector;
  priorUseForBiasVector.resize( PriorNames.size() );
  AtlasDefTable.add(currentRow, 0, "Prior For Bias");
  AtlasDefTable.add(currentRow, 1, ": [");
  AtlasDefTable.add(currentRow, PriorNames.size() + 2 + 1, "]");
  for( unsigned int pwi = 0; pwi < PriorNames.size(); pwi++ )
    {
    priorUseForBiasVector[pwi] = atlasDefinitionParser.GetUseForBias(PriorNames[pwi]);
    AtlasDefTable.add(currentRow, 2 + pwi, priorUseForBiasVector[pwi], "%d");
    }

    { // Print out the ranges.
    currentRow++;
    for( unsigned int pwi = 0; pwi < PriorNames.size(); pwi++ )
      {
      std::map<std::string, AtlasDefinition::BoundsType> temp_range_List;
      for( unsigned int tt = 0; tt < inputVolumeTypes.size(); tt++ )
        {
        AtlasDefTable.add(currentRow + tt * 2 + 0, 0, std::string(inputVolumeTypes[tt]) + std::string(" Lower") );
        AtlasDefTable.add(currentRow + tt * 2 + 0, 1, ": [");
        AtlasDefTable.add(currentRow + tt * 2 + 0, PriorNames.size() + 2 + 1, "]");

        AtlasDefTable.add(currentRow + tt * 2 + 1, 0, std::string(inputVolumeTypes[tt]) + std::string(" Upper") );
        AtlasDefTable.add(currentRow + tt * 2 + 1, 1, ": [");
        AtlasDefTable.add(currentRow + tt * 2 + 1, PriorNames.size() + 2 + 1, "]");

        temp_range_List[inputVolumeTypes[tt]] = atlasDefinitionParser.GetBounds(PriorNames[pwi], inputVolumeTypes[tt]);
        AtlasDefTable.add(currentRow + tt * 2 + 0, 2 + pwi, temp_range_List[inputVolumeTypes[tt]].GetLower(), "%4.2f");
        AtlasDefTable.add(currentRow + tt * 2 + 1, 2 + pwi, temp_range_List[inputVolumeTypes[tt]].GetUpper(), "%4.2f");
        }
      }
    }

    {
    std::ostringstream oss;
    AtlasDefTable.Print(oss);
    muLogMacro( << oss.str() );
    }
  muLogMacro(
    << "Max bias polynomial degree: " << maxBiasDegree << std::endl );
  muLogMacro(<< "Atlas warping: " << !atlasWarpingOff << std::endl );
  muLogMacro(
    << "Atlas warp spline grid size: " << gridSize[0] << " X "
    << gridSize[1] << " X "
    << gridSize[2] << std::endl );
  muLogMacro(<< std::endl );
  muLogMacro(<< "=== Start ===\n");
  muLogMacro(<< "Registering images using affine transform...\n");

  GenericTransformType::Pointer atlasToSubjectPreSegmentationTransform = NULL;

  std::vector<FloatImagePointer> atlasOriginalImageList;
  ByteImagePointer               atlasBrainMask;
    { // Read template images needed for atlas registration
      // muLogMacro(<< "Read template mask");
    const std::string templateMask = atlasDefinitionParser.GetTemplateBrainMask();
    if( templateMask.size() < 1 )
      {
      muLogMacro( <<  "No template mask specified" << std::endl );
      return EXIT_FAILURE;
      }
    typedef itk::ImageFileReader<ByteImageType> ReaderType;
    typedef ReaderType::Pointer                 ReaderPointer;

    muLogMacro( << "Reading mask : " << templateMask << "...\n");

    ReaderPointer imgreader = ReaderType::New();
    imgreader->SetFileName( templateMask.c_str() );

    try
      {
      imgreader->Update();
      }
    catch( ... )
      {
      muLogMacro( << "ERROR:  Could not read image " << templateMask << "." << std::endl );
      return EXIT_FAILURE;
      }
    atlasBrainMask = imgreader->GetOutput();
    }

  std::vector<FloatImagePointer> intraSubjectRegisteredImageList;
  std::vector<FloatImagePointer> intraSubjectRegisteredRawImageList;
  std::vector<std::string>       priorfnlist;
  std::vector<std::string>       templateVolumes( inputVolumeTypes.size() );
  for( unsigned int q = 0; q < inputVolumeTypes.size(); q++ )
    {
    const AtlasDefinition::TemplateMap &         tm = atlasDefinitionParser.GetTemplateVolumes();
    AtlasDefinition::TemplateMap::const_iterator ti = tm.find(inputVolumeTypes[q]);
    if( ti != tm.end() )
      {
      std::string temp = ti->second;
      std::cerr << "STATUS:  Atlas image of type: " << inputVolumeTypes[q] << " added with filename: " << temp
                << std::endl;
      templateVolumes[q] = temp;
      }
    else
      {
      std::cerr << "ERROR:  Atlas image of type: " << inputVolumeTypes[q] << " not found in xml file." << std::endl;
      throw;
      }
    }

  std::vector<bool> duplicatesFound;
    {
    typedef AtlasRegistrationMethod<float, float> AtlasRegType;
    AtlasRegType::Pointer atlasreg = AtlasRegType::New();

    if( debuglevel > 0 )
      {
      atlasreg->DebugOn();
      atlasreg->SetDebugLevel(debuglevel);
      }

    atlasreg->SetSuffix(defaultSuffix);
    // Compute list of file names for the atlasOriginalPriors
    for( unsigned int q = 0; q < PriorNames.size(); q++ )
      {
      priorfnlist.push_back( atlasDefinitionParser.GetPriorFilename( PriorNames[q] ) );
      }
      {
      std::vector<FloatImagePointer> intraSubjectRawImageList;
      intraSubjectRawImageList.clear();
      intraSubjectRawImageList.resize(inputVolumes.size(), 0);
      std::vector<FloatImagePointer> intraSubjectNoiseRemovedImageList;
      intraSubjectNoiseRemovedImageList.clear();
      intraSubjectNoiseRemovedImageList.resize(inputVolumes.size(), 0);
        { // StartOriginalImagesList
        const std::string suffixstr = "";
          { // Read subject images needed for atlas registration
            // muLogMacro(<< "Read subject images");
          if( inputVolumes.size() < 1 )
            {
            muLogMacro( <<  "No data images specified" << std::endl );
            return EXIT_FAILURE;
            }

          typedef itk::ImageFileReader<FloatImageType> ReaderType;
          typedef ReaderType::Pointer                  ReaderPointer;

          std::vector<std::string> intraSubjectTransformFileNames( inputVolumes.size() );
          for( unsigned int i = 0; i < inputVolumes.size(); i++ )
            {
            muLogMacro(
              << "Reading image " << i + 1 << ": " << inputVolumes[i] << "...\n");

            ReaderPointer imgreader = ReaderType::New();
            imgreader->SetFileName( inputVolumes[i].c_str() );

            try
              {
              imgreader->Update();
              }
            catch( ... )
              {
              muLogMacro( << "ERROR:  Could not read image " << inputVolumes[i] << "." << std::endl );
              return EXIT_FAILURE;
              }
            // Initialize with file read in
            FloatImageType::Pointer typewiseEqualizedToFirstImage = imgreader->GetOutput();
#if 0       // This needs more testing.
            // Now go looking to see if this image type has already been found,
            // and equalize to the first image of this type if found.
            for( unsigned int prevImageIndex = 0; prevImageIndex < i; prevImageIndex++ )
              {
              if( inputVolumeTypes[i] == inputVolumeTypes[prevImageIndex] )
              // If it matches a  previous found image type,
              // then histogram equalize
                {
                muLogMacro( << "Equalizing image (" << i << ") to image (" << prevImageIndex << ")" << std::endl );
                typedef itk::HistogramMatchingImageFilter<FloatImageType,
                                                          FloatImageType> HistogramMatchingFilterType;
                HistogramMatchingFilterType::Pointer histogramfilter
                  = HistogramMatchingFilterType::New();

                histogramfilter->SetInput( imgreader->GetOutput() );
                histogramfilter->SetReferenceImage( intraSubjectNoiseRemovedImageList[prevImageIndex] );

                histogramfilter->SetNumberOfHistogramLevels( 128 );
                histogramfilter->SetNumberOfMatchPoints( 16 );
                // histogramfilter->ThresholdAtMeanIntensityOn();
                histogramfilter->Update();
                // Overwrite if necessary.
                typewiseEqualizedToFirstImage = histogramfilter->GetOutput();
                break;
                }
              }
#endif

            // Normalize Image Intensities:
            muLogMacro( << "Standardizing Intensities: ...\n" );
            intraSubjectRawImageList[i] = StandardizeMaskIntensity<FloatImageType, ByteImageType>(
                typewiseEqualizedToFirstImage,
                NULL,
                0.0005, 1.0 - 0.0005,
                1, 0.95 * MAX_IMAGE_OUTPUT_VALUE,
                0, MAX_IMAGE_OUTPUT_VALUE);
            muLogMacro( << "done.\n" );
#if 1
              {
              std::vector<unsigned int> unused_gridSize;
              double                    localFilterTimeStep = filterTimeStep;
              if( localFilterTimeStep <= 0 )
                {
                FloatImageType::SpacingType::ValueType minPixelSize =
                  vcl_numeric_limits<FloatImageType::SpacingType::ValueType>::max();
                const FloatImageType::SpacingType & imageSpacing = intraSubjectRawImageList[i]->GetSpacing();
                for( int is = 0; is < FloatImageType::ImageDimension; ++is )
                  {
                  minPixelSize = vcl_min( minPixelSize, imageSpacing[is]);
                  }
                localFilterTimeStep =
                  ( (minPixelSize - vcl_numeric_limits<FloatImageType::SpacingType::ValueType>::epsilon() )
                    / ( vcl_pow(2.0, FloatImageType::ImageDimension + 1 ) )
                  );
                }
              intraSubjectNoiseRemovedImageList[i] =
                DenoiseFiltering<FloatImageType>(intraSubjectRawImageList[i], filterMethod, filterIteration,
                                                 localFilterTimeStep,
                                                 unused_gridSize);
              if( debuglevel > 1 )
                {
                // DEBUG:  This code is for debugging purposes only;
                typedef itk::ImageFileWriter<FloatImageType> WriterType;
                WriterType::Pointer writer = WriterType::New();
                writer->UseCompressionOn();

                std::stringstream template_index_stream("");
                template_index_stream << i;
                const std::string fn = outputDir + "/DENOISED_INDEX_" + template_index_stream.str() + ".nii.gz";
                writer->SetInput(intraSubjectNoiseRemovedImageList[i]);
                writer->SetFileName(fn.c_str() );
                writer->Update();
                muLogMacro( << "DEBUG:  Wrote image " << fn <<  std::endl);
                }
              }
#else
            intraSubjectNoiseRemovedImageList[i] = intraSubjectRawImageList[i];
#endif
            intraSubjectTransformFileNames[i] = outputDir
              + GetStripedImageFileNameExtension(inputVolumes[i]) + std::string(
                "_to_")
              + GetStripedImageFileNameExtension(inputVolumes[0]) + suffixstr
              + std::string(".mat");
            }
          atlasreg->SetIntraSubjectOriginalImageList(intraSubjectNoiseRemovedImageList);
          atlasreg->SetIntraSubjectTransformFileNames(intraSubjectTransformFileNames);
          }

          { // Read template images needed for atlas registration
            // muLogMacro(<< "Read template images");
          if( templateVolumes.size() < 1 )
            {
            muLogMacro( <<  "No data images specified" << std::endl );
            return EXIT_FAILURE;
            }

          typedef itk::ImageFileReader<FloatImageType> ReaderType;
          typedef ReaderType::Pointer                  ReaderPointer;

          atlasOriginalImageList.clear();
          atlasOriginalImageList.resize(templateVolumes.size(), 0);
          for( unsigned int atlasIndex = 0; atlasIndex < templateVolumes.size(); atlasIndex++ )
            {
            // KENT: HACK:  This currently just checks one previous location in the list to
            //             determine if the image was already loaded,  It should check the
            //             entire list and avoid loading duplicate images into the list
            //             of atlas images.
            if( ( atlasIndex > 0 ) && ( templateVolumes[atlasIndex] == templateVolumes[atlasIndex - 1] ) )
              { // If they are the same name, then just use same reference
              muLogMacro(
                << "Referencing previous image " << atlasIndex + 1 << ": " << templateVolumes[atlasIndex] << "...\n");
              atlasOriginalImageList[atlasIndex] = atlasOriginalImageList[atlasIndex - 1];
              }
            else
              {
              muLogMacro(
                << "Reading image " << atlasIndex + 1 << ": " << templateVolumes[atlasIndex] << "...\n");

              ReaderPointer imgreader = ReaderType::New();
              imgreader->SetFileName( templateVolumes[atlasIndex].c_str() );

              try
                {
                imgreader->Update();
                }
              catch( ... )
                {
                muLogMacro( << "ERROR:  Could not read image " << templateVolumes[atlasIndex] << "." << std::endl );
                return EXIT_FAILURE;
                }

              muLogMacro( << "Standardizing Intensities: ..." );
              FloatImagePointer img_i = StandardizeMaskIntensity<FloatImageType, ByteImageType>(
                  imgreader->GetOutput(),
                  atlasBrainMask,
                  0.0005, 1.0 - 0.0005,
                  1, 0.95 * MAX_IMAGE_OUTPUT_VALUE,
                  0, MAX_IMAGE_OUTPUT_VALUE);
              muLogMacro( << "done." << std::endl );
              atlasOriginalImageList[atlasIndex] = img_i;
              if( debuglevel > 7 )
                {
                typedef itk::ImageFileWriter<FloatImageType> FloatWriterType;
                FloatWriterType::Pointer writer = FloatWriterType::New();

                std::stringstream write_atlas_index_stream("");
                write_atlas_index_stream << atlasIndex;
                const std::string fn
                  = outputDir + std::string("RenormalizedAtlasTemplate_") + write_atlas_index_stream.str() + suffstr;

                writer->SetInput(atlasOriginalImageList[atlasIndex] );
                writer->SetFileName( fn.c_str() );
                writer->UseCompressionOn();
                writer->Update();
                }
              }
            }
          atlasreg->SetAtlasOriginalImageList(atlasOriginalImageList);
          atlasreg->SetInputVolumeTypes(inputVolumeTypes);
          const std::string atlasTransformFileName = outputDir
            + GetStripedImageFileNameExtension(templateVolumes[0])
            + std::string("_to_")
            + GetStripedImageFileNameExtension(inputVolumes[0]) + suffixstr
            + std::string("PreSegmentation.mat");
          atlasreg->SetAtlasToSubjectTransformFileName(atlasTransformFileName);
          }

        //    atlasreg->SetOutputDebugDir(outputDir);

        if( !( ( atlasToSubjectTransformType.compare("Identity") == 0 )
               || ( atlasToSubjectTransformType.compare("Rigid") == 0 )
               || ( atlasToSubjectTransformType.compare("Affine") == 0 )
               || ( atlasToSubjectTransformType.compare("BSpline") == 0 ) )
            )
          {
          muLogMacro(
            "ERROR:  Invalid atlasToSubjectTransformType specified" << atlasToSubjectTransformType << std::endl);
          return EXIT_FAILURE;
          }

        if( !( ( subjectIntermodeTransformType.compare("Identity") == 0 )
               || ( subjectIntermodeTransformType.compare("Rigid") == 0 )
               || ( subjectIntermodeTransformType.compare("Affine") == 0 )
               || ( subjectIntermodeTransformType.compare("BSpline") == 0 ) )
            )
          {
          muLogMacro(
            "ERROR:  Invalid subjectIntermodeTransformType specified" << subjectIntermodeTransformType << std::endl);
          return EXIT_FAILURE;
          }

        if( atlasToSubjectInitialTransform != "" )
          {
          muLogMacro(<< "atlasToSubjectInitialTransform specified." << std::endl)
          if( atlasToSubjectTransformType.compare("Identity") == 0 )
            {
            // Error because we're applying an identity transform by an initial transform was supplied.

            muLogMacro(
              << "ERROR:  atlasToSubjectTransformType is Identity but an initial transform supplied." << std::endl);
            return EXIT_FAILURE;
            }
          try
            {
            GenericTransformType::Pointer atlasToSubjectCurrentGenericTransform = itk::ReadTransformFromDisk(
                atlasToSubjectInitialTransform);

            const std::string initialTransformFileType = atlasToSubjectCurrentGenericTransform->GetNameOfClass();
            try
              {
              if( initialTransformFileType == "VersorRigid3DTransform" )
                {
                if( !( (atlasToSubjectTransformType.compare("Rigid") == 0 )
                       || ( atlasToSubjectTransformType.compare("Affine") == 0 )
                       || ( atlasToSubjectTransformType.compare("BSpline") == 0 ) )
                    )
                  {
                  muLogMacro(
                    << "Error: initialAtlasToSubjectTransform is a VersorRigid3DTransform but atlasToSubhectTransfromType is not Rigid, Affine, or BSpline."
                    << std::endl);
                  return EXIT_FAILURE;
                  }
                }
              else if( initialTransformFileType == "ScaleVersor3DTransform" )
                {
                if( !( (atlasToSubjectTransformType.compare("Rigid") == 0 )
                       || ( atlasToSubjectTransformType.compare("Affine") == 0 )
                       || ( atlasToSubjectTransformType.compare("BSpline") == 0 ) )
                    )
                  {
                  muLogMacro(
                    << "Error: initialAtlasToSubjectTransform is a ScaleVersor3DTransform but atlasToSubhectTransfromType is not Rigid, Affine, or BSpline."
                    << std::endl);
                  return EXIT_FAILURE;
                  }
                }
              else if( initialTransformFileType == "ScaleSkewVersor3DTransform" )
                {
                if( !( (atlasToSubjectTransformType.compare("Rigid") == 0 )
                       || ( atlasToSubjectTransformType.compare("Affine") == 0 )
                       || ( atlasToSubjectTransformType.compare("BSpline") == 0 ) )
                    )
                  {
                  muLogMacro(
                    << "Error: initialAtlasToSubjectTransform is a ScaleSkewVersor3DTransform but atlasToSubhectTransfromType is not Rigid, Affine, or BSpline."
                    << std::endl);
                  return EXIT_FAILURE;
                  }
                }
              else if( initialTransformFileType == "AffineTransform" )
                {
                if( !( ( atlasToSubjectTransformType.compare("Affine") == 0 )
                       || ( atlasToSubjectTransformType.compare("BSpline") == 0 ) )
                    )
                  {
                  muLogMacro(
                    << "Error: initialAtlasToSubjectTransform is a AffineTransform but atlasToSubhectTransfromType is not Affine, or BSpline."
                    << std::endl);
                  return EXIT_FAILURE;
                  }
                }
              else if( initialTransformFileType == "BSplineDeformableTransform" )
                {
                if( !( ( atlasToSubjectTransformType.compare("BSpline") == 0 ) )
                    )
                  {
                  muLogMacro(
                    << "Error: initialAtlasToSubjectTransform is a BSplineDeformableTransform but atlasToSubhectTransfromType is not BSpline."
                    << std::endl);
                  return EXIT_FAILURE;
                  }
                }
              else
                {
                itkGenericExceptionMacro( << "ERROR:  Invalid transform initializer type found:  "
                                          << initialTransformFileType );
                }
              }
            catch( itk::ExceptionObject & excp )
              {
              muLogMacro(<< "Error: error while reading the atlasToSubjectInitialTransform" << std::endl);
              return EXIT_FAILURE;
              }

            atlasreg->SetAtlasToSubjectInitialTransform(atlasToSubjectCurrentGenericTransform);
            }
          catch( itk::ExceptionObject & excp )
            {
            muLogMacro(
              "ERROR:  Invalid atlasToSubjectInitialTransform specified" << atlasToSubjectInitialTransform
                                                                         << std::endl);
            return EXIT_FAILURE;
            }
          }

        atlasreg->SetAtlasLinearTransformChoice(atlasToSubjectTransformType);
        atlasreg->SetImageLinearTransformChoice(subjectIntermodeTransformType);

        atlasreg->SetWarpGrid(
          gridSize[0],
          gridSize[1],
          gridSize[2]);
        muLogMacro(<< "Registering and resampling images..." << std::endl);
        // EMSTimer* regtimer = new EMSTimer();
        //      muLogMacro(<< "Registration took " <<
        // regtimer->GetElapsedHours() << " hours, ");
        //      muLogMacro(<< regtimer->GetElapsedMinutes() << " minutes, ");
        //      muLogMacro(<< regtimer->GetElapsedSeconds() << " seconds\n");
        //      delete regtimer;
        itk::TimeProbe regtimer;
        regtimer.Start();
        atlasreg->SetOutputDebugDir(outputDir);
        try
          {
          atlasreg->Update();
          }
        catch( itk::ExceptionObject & e )
          {
          std::cerr << "Exception caught!" << std::endl;
          std::cerr << e << std::endl;
          return EXIT_FAILURE;
          }

        regtimer.Stop();
        itk::RealTimeClock::TimeStampType elapsedTime
          = regtimer.GetTotal();
        muLogMacro(<< "Registration took " << elapsedTime << " " << regtimer.GetUnit() << std::endl);

        std::vector<GenericTransformType::Pointer> intraSubjectTransforms = atlasreg->GetIntraSubjectTransforms();

        // ::ResampleImages()
          {
          // muLogMacro(<< "ResampleImages");

          // Define the internal reader type
          typedef itk::ResampleImageFilter<FloatImageType, FloatImageType> ResampleType;
          typedef ResampleType::Pointer                                    ResamplePointer;

          intraSubjectRegisteredImageList =
            ResampleImageList(resamplerInterpolatorType, intraSubjectNoiseRemovedImageList,
                              intraSubjectTransforms);
          intraSubjectRegisteredRawImageList = ResampleImageList(resamplerInterpolatorType, intraSubjectRawImageList,
                                                                 intraSubjectTransforms);
          assert( intraSubjectRegisteredImageList.size() == intraSubjectNoiseRemovedImageList.size() );
          assert( intraSubjectRegisteredImageList.size() == intraSubjectRawImageList.size() );
          intraSubjectNoiseRemovedImageList.clear();
          intraSubjectRawImageList.clear();
          } // End registering data

          {
          // Now check that the intraSubjectNoiseRemovedImageList has positive
          // definite covariance matrix.
          // The algorithm is not stable if the covariance matrix is not
          // positive definite, and this
          // occurs when two or more of the images are linearly dependant (i.e.
          // nearly the same image).
          duplicatesFound = FindDuplicateImages(intraSubjectRegisteredImageList);
          for( size_t q = 0; q < duplicatesFound.size(); q++ )
            {
            if( duplicatesFound[q] == true )
              {
              std::cout << "WARNING: Found images that were very highly correlated." << std::endl;
              std::cout << "WARNING: Removing image " << inputVolumes[q] << " from further processing" << std::endl;
              std::cout << "WARNING:" << std::endl;
              }
            }
          }
        } // EndOriginalImagesList

      atlasToSubjectPreSegmentationTransform = atlasreg->GetAtlasToSubjectTransform();
      if( debuglevel > 9 )
        { // NOTE THIS IS REALLY ANNOYING FOR LARGE BSPLINE REGISTRATIONS!
        muLogMacro( << __FILE__ << " " << __LINE__ << " " << atlasToSubjectPreSegmentationTransform->GetFixedParameters(
                      ) << std::endl );
        muLogMacro(
          << __FILE__ << " " << __LINE__ << " " << atlasToSubjectPreSegmentationTransform->GetParameters()
          << std::endl );
        }
      if( debuglevel > 4 )
        {
        // Write the registered template and images
        if( !writeLess )
          {
          muLogMacro(<< "Writing registered template...\n");

          typedef itk::ResampleImageFilter<FloatImageType, FloatImageType> ResampleType;
          typedef ResampleType::Pointer                                    ResamplePointer;
          ResamplePointer resampler = ResampleType::New();

          resampler->SetInput(atlasOriginalImageList[0]);
          resampler->SetTransform(atlasToSubjectPreSegmentationTransform);

          resampler->SetOutputParametersFromImage(intraSubjectRegisteredRawImageList[0]);
          resampler->SetDefaultPixelValue(0);
          resampler->Update();
          typedef itk::CastImageFilter<FloatImageType, ShortImageType>
            ShortRescaleType;
          ShortRescaleType::Pointer rescaler = ShortRescaleType::New();
          rescaler->SetInput( resampler->GetOutput() );
          rescaler->Update();

          typedef itk::ImageFileWriter<ShortImageType> ShortWriterType;
          ShortWriterType::Pointer writer = ShortWriterType::New();

          std::string fn
            = outputDir + std::string("AtlasToSubjectInitialization_") + suffstr;

          writer->SetInput( rescaler->GetOutput() );
          writer->SetFileName( fn.c_str() );
          writer->UseCompressionOn();
          writer->Update();
          }
        for( unsigned int i = 0; i < intraSubjectRegisteredRawImageList.size(); i++ )
          {
          typedef itk::RescaleIntensityImageFilter<FloatImageType, ShortImageType>
            ShortRescaleType;

          ShortRescaleType::Pointer rescaler = ShortRescaleType::New();
          rescaler->SetOutputMinimum(0);
          rescaler->SetOutputMaximum(MAX_IMAGE_OUTPUT_VALUE);
          rescaler->SetInput(intraSubjectRegisteredRawImageList[i]);
          rescaler->Update();

          std::string fn
            = outputDir + GetStripedImageFileNameExtension( ( inputVolumes[i] ) )
              + std::string("_registered") + suffstr;

          typedef itk::ImageFileWriter<ShortImageType> ShortWriterType;
          ShortWriterType::Pointer writer = ShortWriterType::New();

          writer->SetInput( rescaler->GetOutput() );
          writer->SetFileName( fn.c_str() );
          writer->UseCompressionOn();
          writer->Update();
          }
        }
      }
    } // end atlas reg block

    {
    // Now based on the duplicates found, update input lists to only use unique
    // images
    intraSubjectRegisteredImageList = RemoveDuplicates(intraSubjectRegisteredImageList, duplicatesFound);
    intraSubjectRegisteredRawImageList = RemoveDuplicates(intraSubjectRegisteredRawImageList, duplicatesFound);
    inputVolumes = RemoveDuplicates(inputVolumes, duplicatesFound);
    inputVolumeTypes = RemoveDuplicates(inputVolumeTypes, duplicatesFound);
    }

  muLogMacro(<< "Rescale intensity of filtered images...\n");
    {
    RescaleFunctionLocal(intraSubjectRegisteredImageList);
    RescaleFunctionLocal(intraSubjectRegisteredRawImageList);
    }
  // Start the segmentation process.
  for( unsigned int segmentationLevel = 0; segmentationLevel < 1; segmentationLevel++ )
    {
#if 0
    // Prior Names
    // Prior Names         : [ AIR  GM    BGM  CSF  NOTCSF NOTGM NOTVB NOTWM VB
    //   WM    ]
    // Prior weight scales : [ 1.00 1.00  1.00 1.00 1.00   1.00  1.00  1.00
    //  1.00 0.50  ]
    // Prior Label Codes   : [ 0    2     3    4    6      7     9     8     5
    //    1     ]
    std::map<unsigned int, std::vector<unsigned int> > ParentLabels;
    // Air
    ParentLabels[100].push_back(0);
    // GM-BGM-NOTGM
    ParentLabels[101].push_back(2);
    ParentLabels[101].push_back(3);
    ParentLabels[101].push_back(7);
    // CSF-NOTCSF
    ParentLabels[102].push_back(4);
    ParentLabels[102].push_back(6);
    // WM-NOTWM
    ParentLabels[103].push_back(1);
    ParentLabels[103].push_back(8);
    // VB-NOTVB
    ParentLabels[104].push_back(5);
    ParentLabels[104].push_back(9);

    std::map<unsigned int, std::string> ParentLabelNames;
    for(
      std::map<unsigned int, std::vector<unsigned int> >::const_iterator plIter = ParentLabels.begin();
      plIter != ParentLabels.end(); plIter++ )
      {
      std::string combinedName = "";
      bool        firstpass = true;
      for(
        std::vector<unsigned int>::const_iterator vIter = plIter->second.begin();
        vIter != plIter->second.end(); vIter++ )
        {
        if( firstpass )
          {
          firstpass = false;
          }
        else
          {
          combinedName += "-";
          }
        combinedName += PriorNames[*vIter];
        }
      ParentLabelNames[plIter->first] = combinedName;
      }

#endif

    SegFilterType::Pointer segfilter = SegFilterType::New();
    // __MAX__PROBS
      {
      std::vector<FloatImagePointer> atlasOriginalPriors( PriorNames.size() );
      for( unsigned int i = 0; i < PriorNames.size(); i++ )
        {
        typedef itk::ImageFileReader<FloatImageType> ReaderType;
        typedef ReaderType::Pointer                  ReaderPointer;
        ReaderPointer priorReader = ReaderType::New();
        priorReader->SetFileName( atlasDefinitionParser.GetPriorFilename(PriorNames[i]) );
        priorReader->Update();
        FloatImageType::Pointer temp = priorReader->GetOutput();
        atlasOriginalPriors[i] = temp;
        }
      segfilter->SetPriors(atlasOriginalPriors);
      }
      {
      SegFilterType::RangeDBType myRanges;
      for( unsigned int q = 0; q < PriorNames.size(); q++ )
        {
        std::map<std::string, AtlasDefinition::BoundsType> temp_range_List;
        for( unsigned int tt = 0; tt < inputVolumeTypes.size(); tt++ )
          {
          temp_range_List[inputVolumeTypes[tt]] = atlasDefinitionParser.GetBounds(PriorNames[q], inputVolumeTypes[tt]);
          }
        myRanges[PriorNames[q]] = temp_range_List;
        }
      segfilter->SetTissueTypeThresholdMapsRange(myRanges);
      segfilter->SetPriorNames(PriorNames);
      }

    // //Static component that does not depend on priors-consolidation.
    muLogMacro(<< "Start segmentation...\n");
    segfilter->SetOutputDebugDir(outputDir);

    if( debuglevel > 0 )
      {
      segfilter->DebugOn();
      segfilter->SetDebugLevel(debuglevel);
      }

    segfilter->SetInputImages(intraSubjectRegisteredImageList);
    segfilter->SetRawInputImages(intraSubjectRegisteredRawImageList);
    segfilter->SetInputVolumeTypes(inputVolumeTypes);

    segfilter->SetMaximumIterations( maxIterations );
    segfilter->SetOriginalAtlasImages(atlasOriginalImageList);
    segfilter->SetTemplateBrainMask(atlasBrainMask);
    segfilter->SetTemplateGenericTransform(atlasToSubjectPreSegmentationTransform);

    segfilter->SetPriorWeights(priorsWeightList);
    segfilter->SetPriorLabelCodeVector(priorLabelCodeVector);
    segfilter->SetPriorUseForBiasVector(priorUseForBiasVector);
    segfilter->SetPriorIsForegroundPriorVector(priorIsForegroundPriorVector);

    segfilter->SetMaxBiasDegree(maxBiasDegree);

#if 0
    // THIS ACTAULLY NEEDS TO BE USED FOR ONLY DOING THE FINAL ITERATION OF THE
    // SEGMENTATION PROCESS AND SKIPPING MOST OF EMLoop
    // If this "Post" transform is given, then jump over the looping, and
    // warp priors,do estimates, and bias correct in one step.
      {
      ReadGenericTransform;

      segfilter->SetTemplateGenericTransform();
      // Turn off warping, and just use the input given from disk.
      }
#endif

    if( !atlasWarpingOff )
      {
      segfilter->UpdateTransformationOn();
      }
    else
      {
      segfilter->UpdateTransformationOff();
      }
    segfilter->SetWarpGrid(
      gridSize[0],
      gridSize[1],
      gridSize[2]);

    segfilter->Update();

    std::vector<std::string> names = inputVolumes;
    // Write the secondary outputs
    if( !writeLess )
      {
      muLogMacro(<< "Writing filtered and bias corrected images...\n");
      // std::vector<FloatImagePointer> imgset = segfilter->GetCorrected();
      std::vector<FloatImagePointer> imgset = segfilter->GetRawCorrected();
      std::vector<std::string>       outFileNames;
      if( outputVolumes.size() == 1 )
        {
        outFileNames.resize(imgset.size() );
        for( unsigned i = 0; i < imgset.size(); i++ )
          {
          char buf[8192];
          sprintf(buf, outputVolumes[0].c_str(), inputVolumeTypes[i].c_str(), i);
          outFileNames[i] = buf;
          }
        }
      else if( outputVolumes.size() != imgset.size() )
        {
        std::cerr << imgset.size() << " images in filter output, but "
                  << outputVolumes.size() << " names in output volumes list"
                  << std::endl;
        return EXIT_FAILURE;
        }
      else
        {
        outFileNames = outputVolumes;
        }
      for( unsigned i = 0; i < imgset.size(); i++ )
        {
        // typedef itk::RescaleIntensityImageFilter<FloatImageType,
        // ShortImageType> ShortRescaleType;
        typedef itk::CastImageFilter<FloatImageType, ShortImageType> CasterType;
        CasterType::Pointer caster = CasterType::New();

        caster->SetInput(imgset[i]);
        caster->Update();
        // std::string fn
        //  = outputDir + GetStripedImageFileNameExtension(names[i]) + std::string("_corrected")
        //    + suffstr;

        typedef itk::ImageFileWriter<ShortImageType> ShortWriterType;
        ShortWriterType::Pointer writer = ShortWriterType::New();

        writer->SetInput( caster->GetOutput() );
        writer->SetFileName(outFileNames[i]);
        writer->UseCompressionOn();
        writer->Update();
        }
      }

    // Average together all the input images of a given type
    typedef std::map<std::string, std::vector<FloatImageType::Pointer> > imgTypeMapType;
    imgTypeMapType volumesByImgType;

    std::vector<FloatImagePointer> imgset = segfilter->GetRawCorrected();
    for( unsigned int k = 0; k < inputVolumeTypes.size(); ++k )
      {
      volumesByImgType[inputVolumeTypes[k]].push_back(imgset[k]);
      }

    imgTypeMapType::const_iterator mapItr;
    for( mapItr = volumesByImgType.begin();
         mapItr != volumesByImgType.end();
         ++mapItr )
      {
      std::string volumeType = mapItr->first;
      // Can't average images of type other since it's really a mix of types.
      std::transform(volumeType.begin(), volumeType.end(), volumeType.begin(), tolower);
      if( volumeType == "other" )
        {
        continue;
        }

      FloatImageType::Pointer avgImage = AverageImageList(mapItr->second);
      // Write out average image.
      typedef itk::CastImageFilter<FloatImageType, ShortImageType> CasterType;
      CasterType::Pointer caster = CasterType::New();

      caster->SetInput(avgImage);
      caster->Update();

      std::string avgFileName = outputDir + volumeType + std::string("_average") + suffstr;

      muLogMacro(<< "Writing averaged corrected input images... " << avgFileName << std::endl );

      typedef itk::ImageFileWriter<ShortImageType> ShortWriterType;
      ShortWriterType::Pointer writer = ShortWriterType::New();

      writer->SetInput( caster->GetOutput() );
      writer->SetFileName(avgFileName);
      writer->UseCompressionOn();
      writer->Update();
      }

    // Write warped template and bspline trafo
    if( !atlasWarpingOff )
      {
      std::vector<SegFilterType::InputImagePointer> WarpedAtlasList = segfilter->GenerateWarpedAtlasImages();
      for( unsigned int index = 0; index < WarpedAtlasList.size(); index++ )
        {
        typedef itk::RescaleIntensityImageFilter<FloatImageType, ByteImageType>
          ByteRescaleType;

        ByteRescaleType::Pointer rescaler = ByteRescaleType::New();
        rescaler->SetOutputMinimum(0);
        rescaler->SetOutputMaximum(255);
        rescaler->SetInput(WarpedAtlasList[index]);
        rescaler->Update();

        typedef itk::ImageFileWriter<ByteImageType> ByteWriterType;
        ByteWriterType::Pointer writer = ByteWriterType::New();

        const std::string fn = outputDir
          + GetStripedImageFileNameExtension(templateVolumes[index]) + std::string("_to_")
          + GetStripedImageFileNameExtension( ( inputVolumes[0] ) )
          + std::string("_warped") + std::string(".nii.gz");

        muLogMacro(<< "Writing warped template images... " << fn << std::endl );
        writer->SetInput( rescaler->GetOutput() );
        writer->SetFileName( fn.c_str() );
        writer->UseCompressionOn();
        writer->Update();
        }
      if( atlasToSubjectTransform != "" )
        {
        /* HACK Remove this completely
      const std::string postSegmentationTransformFileName = outputDir
        + GetStripedImageFileNameExtension(templateVolumes[0])
        + std::string("_to_")
        + GetStripedImageFileNameExtension( ( inputVolumes[0] ) )
        + std::string("_") + defaultSuffix + "_PostSegmentation.mat";
        */
        const std::string postSegmentationTransformFileName = atlasToSubjectTransform;
        // NOTE:  Aliasing of smart-pointers up the polymorphic tree OK here
        // because
        // the primary
        // smart pointer is gaunteed to exists during lifespan of all aliases.
        muLogMacro(<< "Writing final atlas to subject template... " << postSegmentationTransformFileName << std::endl );
        GenericTransformType::Pointer atlasToSubjectPostSegmentationTransform =
          segfilter->GetTemplateGenericTransform();
        WriteTransformToDisk(atlasToSubjectPostSegmentationTransform, postSegmentationTransformFileName);
        // TODO:  Need to write a short circuit so that if this final transform
        // filename exists, it will just read it in and use it directly
        // without doing all the iterations.
        }
      }

    // Write the labels
    muLogMacro(<< "Writing labels...\n");
      {
      typedef itk::ImageFileWriter<ByteImageType> ByteWriterType;
      ByteWriterType::Pointer writer = ByteWriterType::New();

      writer->SetInput( segfilter->GetOutput() );
      std::string fn;
      if( outputLabels == "" )
        {
        fn = outputDir;
        fn += GetStripedImageFileNameExtension(names[0]);
        fn += std::string("_DirtyLabels");
        fn += suffstr;
        }
      else
        {
        fn = outputDirtyLabels;
        }
      writer->SetFileName( fn.c_str() );
      writer->UseCompressionOn();
      writer->Update();
      }
      {
      typedef itk::ImageFileWriter<ByteImageType> ByteWriterType;
      ByteWriterType::Pointer writer = ByteWriterType::New();

      std::string fn;
      if( outputLabels == "" )
        {
        fn = outputDir;
        fn += GetStripedImageFileNameExtension(names[0]);
        fn += std::string("_labels");
        fn += suffstr;
        }
      else
        {
        fn = outputLabels;
        }

      writer->SetInput( segfilter->GetCleanedOutput() );
      writer->SetFileName( fn.c_str() );
      writer->UseCompressionOn();

      try
        {
        writer->Update();

        fn = outputDir;
        fn += "thresholded_labels.nii.gz";
        writer->SetInput( segfilter->GetThresholdedOutput() );
        writer->SetFileName( fn );
        writer->UseCompressionOn();
        writer->Modified();
        writer->Update();
        }
      catch( itk::ExceptionObject & e )
        {
        muLogMacro( <<  e << std::endl );
        return -1;
        }
      catch( std::exception & e )
        {
        muLogMacro( <<  "Exception: " << e.what() << std::endl );
        return -1;
        }
      catch( std::string & s )
        {
        muLogMacro( <<  "Exception: " << s << std::endl );
        return -1;
        }
      catch( ... )
        {
        muLogMacro( <<  "Unknown exception" << std::endl );
        muLogMacro( << "failed to write image " << fn << std::endl );
        return -1;
        }
      }
    // Write final Posteriors
      {
      // NOTE :  Priors and Posteriors should correspond, so use the PriorNames
      // to get the names.
      for( unsigned int probabilityIndex = 0; probabilityIndex < PriorNames.size(); probabilityIndex++ )
        {
        std::string fn;
        if( posteriorTemplate == "" )
          {
          fn = outputDir;
          fn += GetStripedImageFileNameExtension(names[0]);
          fn += "_POSTERIOR_";
          fn += PriorNames[probabilityIndex];
          fn += suffstr;
          }
        else
          {
          char buf[8192];
          sprintf(buf, posteriorTemplate.c_str(),
                  PriorNames[probabilityIndex].c_str() );
          fn = buf;
          }
        typedef itk::ImageFileWriter<FloatImageType> FloatWriterType;
        FloatWriterType::Pointer writer = FloatWriterType::New();

        FloatImageType::Pointer currPosterior = segfilter->GetPosteriors()[probabilityIndex];
        writer->SetInput( currPosterior );
        writer->SetFileName( fn.c_str() );
        writer->UseCompressionOn();
        writer->Update();
        }
      }
    }
  timer.Stop();
  muLogMacro(<< "All segmentation processes took " << timer.GetTotal() << " " << timer.GetUnit() << std::endl );
  return EXIT_SUCCESS;
}
