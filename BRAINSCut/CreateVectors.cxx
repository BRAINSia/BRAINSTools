#include <ProcessDescription.h>
#include <NeuralParams.h>
#include <ANNParams.h>
#include <SVMParams.h>
#include <itksys/SystemTools.hxx>
#include <list>
#include <vector>
#include <iostream>
#include <itkIO.h>
#include "Utilities.h"
#include "XMLParser.h"
#include "Parser.h"
#include "itkImageRandomConstIteratorWithIndex.h"
#include "itkHistogramMatchingImageFilter.h" // KEY_HE
#include "ShuffleVectors.h"

// / Generate Input Output Vectors
/** Training vectors are generated for one image.
 * \ingroup Step2
 * \param CurrASDef AS Description
 * \param ProbabilityMapName Probability Map Filename
 * \param StructureID Structure ID
 * \param ParamName Registration Parameter Filename
 * \param AtlasType Atlas Image Type (T1, T2, ..)
 * \param AtlasFilename Atlas Image Filename
 * \param AtlasLandmark Atlas Landmark Filename
 * \param DefGaussian Image Smoothing Value (applied to images prior to registration)
 * \param ImageTypeList List of Image Types (T1, T2, ..)
 * \param AtlasToSubjRegistrationFilename Current Deformation Field Filename
 * \param InputVectorSize Size of Input Vectors
 * \param GaussianValue Image Smoothing Value (for Registration)
 * \param GradientProfileSize Size of Iris Information in Input Vector
 * \param IrisSize Size of Iris Information in Input Vector
 * \param GaussianValue Image Smoothing Value (for Registration)
 * \param doANN Boolean to Collect ANN Training Vectors
 * \param ANNTrainingSize ANN Training Vector Size
 * \param ANNInputVectorList ANN Input Training Vector List
 * \param ANNOutputVectorList ANN Output Training Vector List
 * \param doSVM Boolean to Collect SVM Training Vectors
 * \param SVMTrainingSize SVM Training Vector Size
 * \param SVMInputVectorList SVM Input Training Vector List
 * \param SVMOutputVectorList SVM Output Training Vector List
 */
static bool GenerateInputOutputVectors
(
  DataSet *curDataSet,
  const std::string & CurrProbabilityMapName,
  const std::string & RhoMapName,
  const std::string & PhiMapName,
  const std::string & ThetaMapName,
  const std::string & StructureID,
  ProbabilityMapImageType::Pointer DeformedProbabilityMap[],
  DataSet * atlasDataSet, // KEY_HE
  DataSet::TypeVector & ImageTypeList,
  const std::string & AtlasToSubjRegistrationFilename,
  const int InputVectorSize,
  const int GradientProfileSize,
  const float GaussianValue,
  bool doANN,
  unsigned int & ANNVectorsCreated,
  std::ofstream & ANNVectorStream,
  bool doSVM,
  const unsigned int SVMTrainingSize,
  unsigned int & SVMVectorsCreated,
  std::ofstream & SVMVectorStream,
  int currStructureIndex,
  const int NumberOfProbabilityMaps,
  bool histogramEqualization,
  int verbose)
{
  if( !itksys::SystemTools::FileExists( CurrProbabilityMapName.c_str() ) )
    {
    std::cout << "Missing Probability Map File  " << CurrProbabilityMapName
              << std::endl;
    exit(-1);
    }

  const std::string ImageID( curDataSet->GetAttribute<StringValue>("Name") );
  const std::string KnownOutputMaskName( curDataSet->GetMaskFilenameByType(
                                           StructureID) );
  if( KnownOutputMaskName == "" )
    {
    std::cout << "Current KnownOutputMask cannot be found, so Output Vectors Cannot be computed."
              << "\n Skip This Data Set's Adding Vector for Mask :: "
              << StructureID << " " <<     ImageID <<      std::endl;
    return false;
    }
  const std::string landmarkType("");

  // Reading Subject Images....(T1 T2...)
  std::map<std::string, ProbabilityMapImageType::Pointer> MapOfSubject;
  std::map<std::string, ProbabilityMapImageType::Pointer> MapOfAtlas;   // KEY_HE
  std::map<std::string,
           ImageLinearInterpolatorType::Pointer> MapOfImageInterpolators;
  for( DataSet::TypeVector::const_iterator ImageList = ImageTypeList.begin();
       ImageList != ImageTypeList.end();
       ++ImageList )
    {
    const std::string
      curFilename( curDataSet->GetImageFilenameByType(*ImageList) );

    // HACK:  radius is hardcoded to size 0, it should be taken from command
    // line arguments.
    ProbabilityMapImageType::SizeType radius;
    radius[0] = 0; radius[1] = 0; radius[2] = 0;
    MapOfSubject[*ImageList] =
      ReadMedianFilteredImage<ProbabilityMapImageType>(curFilename, radius);
    ImageLinearInterpolatorType::Pointer CurrInterpolator =
      ImageLinearInterpolatorType::New();
    CurrInterpolator->SetInputImage(MapOfSubject[*ImageList]);
    MapOfImageInterpolators[*ImageList] = CurrInterpolator;
    if( verbose > 0 )
      {
      std::cout << " * Current File name :: "
                << curFilename
                << " * Image Origin is :: "
                <<  MapOfSubject[*ImageList]->GetOrigin()
                << std::endl;
      }
    // * Histogram Equalization [ HE ] from Subject to Atlas. KEY_HE
    if( histogramEqualization )
      {
      std::cout << "- Histogram Equalization Performs on the Images " << *ImageList
                << std::endl;
      std::string atlasVolumeFileName = atlasDataSet->GetImageFilenameByType( *ImageList);
      // Check if altlas image is given for the type
      std::transform( atlasVolumeFileName.begin(),
                      atlasVolumeFileName.end(),
                      atlasVolumeFileName.begin(), ::tolower);

      if( atlasVolumeFileName == "na" )
        {
        std::cout << "  !!!!!!!!!!!! CAUTION !!!!!!!!!!!!!!!!!!! " << std::endl
                  << "  -  If atlas data set of the " << *ImageList << std::endl
                  << "     is not given, Histogram Equalization Process will be skipped " << std::endl
                  << "  !!!!!!!!!!!! CAUTION !!!!!!!!!!!!!!!!!!! " << std::endl;
        break;
        }
      MapOfAtlas[*ImageList] =  ReadMedianFilteredImage<ProbabilityMapImageType>(
          atlasDataSet->GetImageFilenameByType( *ImageList), radius );
      typedef itk::HistogramMatchingImageFilter<ProbabilityMapImageType,
                                                ProbabilityMapImageType> HEFilterType;
      HEFilterType::Pointer HEFilter = HEFilterType::New();
      HEFilter->SetReferenceImage( MapOfAtlas[*ImageList] );
      HEFilter->SetInput( MapOfSubject[*ImageList] );
      HEFilter->SetNumberOfHistogramLevels( 100);
      HEFilter->SetNumberOfMatchPoints( 15);
      HEFilter->ThresholdAtMeanIntensityOn();   // TOOD:: Check if we want to do
                                                // this.

      HEFilter->Update();

      // overwrite the input image to the Histogram Equalized one.
      MapOfSubject[*ImageList] = HEFilter->GetOutput();
      if( verbose > 5 )
        {
        itkUtil::WriteImage<ProbabilityMapImageType>(
          MapOfSubject[*ImageList],
          curFilename + "HE.nii.gz" );
        }
      }
    }
  // Use first image in list to get spacing and size
  ProbabilityMapImageType::Pointer           ReferenceImage = MapOfSubject[*( ImageTypeList.begin() )];
  const ProbabilityMapImageType::SpacingType ImageSpacing =
    ReferenceImage->GetSpacing();
  const ProbabilityMapImageType::PointType ImageOrigin =
    ReferenceImage->GetOrigin();
  const ProbabilityMapImageType::SizeType ImageSize =
    ReferenceImage->GetLargestPossibleRegion().
    GetSize();
  itk::ContinuousIndex<float, 3> ImageOuterBoundaryLocationInMM;
  for( int i = 0; i < 3; i++ )
    {
    ImageOuterBoundaryLocationInMM[i] = ImageSize[i] * ImageSpacing[i]
      + ImageOrigin[i];
    }

  if( DeformedProbabilityMap[currStructureIndex].IsNull() )
    {
    std::cout << "No transformation " << AtlasToSubjRegistrationFilename
              << " filename found. Cannot compute registration." << std::endl;
    exit(-1);
    }
  if( verbose > 8 )
    {
    std::cout << "Writing DeformedProbabilityMap with currStructureIndex of "
              << currStructureIndex << std::endl;
    std::cout << "Image Origin is :: "
              << DeformedProbabilityMap[currStructureIndex]->GetOrigin()
              << std::endl;

    itkUtil::WriteImage<RealImageType>(
      DeformedProbabilityMap[currStructureIndex],
      "./debug_def_probmap.nii.gz");
    }
  VerifyNonZeroImage<RealImageType>(DeformedProbabilityMap[currStructureIndex],
                                    CurrProbabilityMapName);

  itk::GradientImageFilter<ProbabilityMapImageType, float,
                           float>::Pointer GradientFilter =
    itk::GradientImageFilter<ProbabilityMapImageType, float, float>::New();
  /*---REGINA:: check if we need this Gradient Filter
                Remove Rho/Phi/Theta and
                  :replace those with parametric variable || get name from xml document
  ---*/
  GradientFilter->SetInput(DeformedProbabilityMap[currStructureIndex]);
  GradientFilter->Update();
  itk::Image<itk::CovariantVector<float,
                                  3>, 3>::Pointer ProbMapGradient = GradientFilter->GetOutput();
  // REGINA:: Read Rho, phi, theta file name from xml
  // Check if rho,phi and theta file exists.
  RealImageType::Pointer DeformedRhoMap =  ImageWarper<RealImageType>(
      AtlasToSubjRegistrationFilename,
      RhoMapName, ReferenceImage);
  RealImageType::Pointer DeformedPhiMap =  ImageWarper<RealImageType>(
      AtlasToSubjRegistrationFilename,
      PhiMapName, ReferenceImage);
  RealImageType::Pointer DeformedThetaMap = ImageWarper<RealImageType>(
      AtlasToSubjRegistrationFilename,
      ThetaMapName, ReferenceImage);
  if( verbose > 0 )
    {
    itkUtil::WriteImage<RealImageType>(
      DeformedRhoMap,
      AtlasToSubjRegistrationFilename
      + "def_rho.nii.gz");
    itkUtil::WriteImage<RealImageType>(
      DeformedPhiMap,
      AtlasToSubjRegistrationFilename
      + "def_phi.nii.gz");
    itkUtil::WriteImage<RealImageType>(
      DeformedThetaMap,
      AtlasToSubjRegistrationFilename
      + "def_theta.nii.gz");
    }

  itk::Index<3> min, max;
  DefineBoundingBox(DeformedProbabilityMap[currStructureIndex], min, max);
  if( min[0] > max[0] || min[1] > max[1] || min[2] > max[2] )
    {
    std::cout
      << "ERROR IN BOUNDING BOX!  Min can not be greater than max!!!"
      << __FILE__ << " " << __LINE__ << std::endl;
    exit(-1);
    }
  ProbabilityMapImageType::RegionType BBregion;
  itk::Size<3>                        BBRange;
  itk::Index<3>                       BBIndex;
  unsigned long                       range[3] = { max[0] - min[0], max[1] - min[1], max[2] - min[2]};
  long int                            index[3] = { min[0], min[1], min[2]};
  BBRange.SetSize(range);
  BBregion.SetSize(BBRange);
  BBIndex.SetIndex(index);
  BBregion.SetIndex(BBIndex);
  itk::ImageRegionIterator<ProbabilityMapImageType> bbri(DeformedProbabilityMap[currStructureIndex],
                                                         BBregion);
  std::vector<float> inputvector(InputVectorSize);
  std::vector<float> outputvector(NumberOfProbabilityMaps);

  unsigned int NonZeroProbMapPoints = 0;
    {
    // Raster throught all voxels within the bounding box of the deformed
    // probability map.
    while( !bbri.IsAtEnd() )
      {
      const ProbabilityMapImageType::IndexType CurrentIndex = bbri.GetIndex();
      ProbabilityMapImageType::PixelType       current_value =
        DeformedProbabilityMap[currStructureIndex]->GetPixel(CurrentIndex);
      if( ( current_value <
            ( HUNDRED_PERCENT_VALUE - PERCENT_MIN_MAX_TOLERANCE ) )
          &&  ( current_value > ( 0.0F + PERCENT_MIN_MAX_TOLERANCE ) )
          )
        {
        NonZeroProbMapPoints++;
        }
      ++bbri;
      }

    std::cout << "Found " << NonZeroProbMapPoints
              << " possible voxels in probability map." << std::endl;
    if( NonZeroProbMapPoints == 0 )
      {
      std::cout << "ERROR:  Could not find any valid probable points for "
                << AtlasToSubjRegistrationFilename + "def_probmap.nii.gz"
                << std::endl;
      exit(-1);
      }
    }
  ProbabilityMapImageType::Pointer KnownOutputMask =
    itkUtil::ReadImage<ProbabilityMapImageType>(KnownOutputMaskName);
  KnownOutputMask =
    itkUtil::ScaleAndCast<ProbabilityMapImageType, ProbabilityMapImageType>(
      KnownOutputMask,
      0,
      HUNDRED_PERCENT_VALUE);
  std::cout
    << "&&&&&&&&&&&&&&&&&&&&&&&&& Smoothing KnownOutputMask by GaussianValue"
    << KnownOutputMaskName << " " << GaussianValue << std::endl;
  if( GaussianValue != 0 )
    {
    KnownOutputMask = SmoothSingleImage<ProbabilityMapImageType>(
        KnownOutputMask,
        GaussianValue);
    itkUtil::WriteImage<ProbabilityMapImageType>(
      KnownOutputMask,
      KnownOutputMaskName
      + "SmoothedMask.nii.gz");
    }

  if( doANN == true )
    {
    // reducer is needed so that we don't run out of memory in creating the
    // maps.
    // TODO:  The MaxNumberOfSamplingPoints per image should be set in the XML
    // file.
    const unsigned int HardMaxNumberSamplingPoints = 500000; // limit to 1500000
                                                             // sampling points
                                                             // per image.
    ProbabilityMapImageType::Pointer ImprobableLocations =
      ProbabilityMapImageType::New();

    CreateNewImageFromTemplate(ImprobableLocations, KnownOutputMask);

    itk::ImageRandomConstIteratorWithIndex<ProbabilityMapImageType>
    randomWalkRegionIterator(DeformedProbabilityMap[currStructureIndex], BBregion);

    randomWalkRegionIterator.SetNumberOfSamples( BBregion.GetNumberOfPixels() );
    // SetNumber of samples to entire region, other methods will cause loop to
    // break before this is hit.
    randomWalkRegionIterator.GoToBegin();
      {
      int SamplesLeftToFind =
        ( HardMaxNumberSamplingPoints <
          NonZeroProbMapPoints ) ?      HardMaxNumberSamplingPoints :
        NonZeroProbMapPoints;
      // Need big buffer to minimize the number of writes in order to speed up
      // the generation of vectors.
      const unsigned int vector_size = outputvector.size()
        + inputvector.size() + 1;
      float *write_buffer =
        new float[SamplesLeftToFind * vector_size];
      unsigned int currentBufferCounter = 0;
      unsigned int value_index = 0;
      while( ( !randomWalkRegionIterator.IsAtEnd() ) && ( SamplesLeftToFind ) )
        {
        const ProbabilityMapImageType::IndexType CurrentIndex =
          randomWalkRegionIterator.GetIndex();
        ProbabilityMapImageType::PixelType current_value =
          static_cast<ProbabilityMapImageType::PixelType>(
            DeformedProbabilityMap[currStructureIndex]
            ->GetPixel(
              CurrentIndex) );
        // cppcheck-suppress variableScope
        neural_scalar_type guard(AUTOSEG_VEC_SENTINEL);
        if( ( current_value <
              ( HUNDRED_PERCENT_VALUE - PERCENT_MIN_MAX_TOLERANCE ) )
            &&  ( current_value > ( 0.0F + PERCENT_MIN_MAX_TOLERANCE ) ) )
          {
          // Generate and add Input Vector
          AddInputVector(inputvector,
                         //  DeformedProbabilityMap[currStructureIndex],
                         ProbMapGradient,
                         DeformedRhoMap,
                         DeformedPhiMap,
                         DeformedThetaMap,
                         //  ImageTypeList,
                         MapOfSubject,
                         MapOfImageInterpolators,
                         DeformedProbabilityMap,
                         NumberOfProbabilityMaps,
                         CurrentIndex,
                         GradientProfileSize /*,                 IrisSize  //We
                                               are not using Iris Size any
                                               More*/
                         );
          // Output Vector
          int structureTrack = 0;
          for( std::vector<float>::iterator fi = outputvector.begin();
               fi != outputvector.end();
               ++fi )
            {
            if( structureTrack == currStructureIndex )
              {
              *fi = static_cast<float>( KnownOutputMask->GetPixel(CurrentIndex) );
              // Ensure that taking values either zero or one
              if( *fi > 0.5F )
                {
                *fi = 1.0F;
                }
              else
                {
                *fi = 0.0F;
                }
              }
            else
              {
              *fi = 0.0F;
              }
            structureTrack++;
            }
          // std::cerr << "output vector size " << outputvector.size() << "
          // input vector size " << inputvector.size() << std::endl;
          // Write out to the file.
          // std::cout << "write output buffer value ::";
          for( std::vector<float>::const_iterator ofi = outputvector.begin();
               ofi != outputvector.end();
               ++ofi )
            {
            // ANNVectorStream.write((const char *)(&(*ofi)), sizeof(*ofi));
            // std::cout  << *ofi<<"   ";
            write_buffer[value_index] = *ofi;
            value_index++;
            }
          // std::cout<<std::endl;
          // std::cout << "write input buffer value ::";
          for( std::vector<float>::const_iterator ifi = inputvector.begin();
               ifi != inputvector.end();
               ++ifi )
            {
            // std::cout <<*ifi<<"   ";
            // ANNVectorStream.write((const char *)(&(*ifi)), sizeof(*ifi));
            write_buffer[value_index] = *ifi;
            value_index++;
            }
          // std::cout<<std::endl;
            { /*Make an Error Image of improbable values*/
            const float ProbError = vcl_abs(
                outputvector[currStructureIndex] - inputvector[0]);
            if( ProbError > 0.66 )
              {
              ImprobableLocations->SetPixel( CurrentIndex,
                                             static_cast<ProbabilityMapImageType::PixelType>( ProbError
                                                                                              *
                                                                                              HUNDRED_PERCENT_VALUE ) );
              }
            }
          // write out sentinel value
          // ANNVectorStream.write((const char *)(&guard), sizeof(guard));
          write_buffer[value_index] = guard;
          // std::cout<<"write guard:: "<<write_buffer[value_index]<<"\n";
          value_index++;
          currentBufferCounter++;
          SamplesLeftToFind--;
          }
        else
          {
          ImprobableLocations->SetPixel(CurrentIndex, current_value);
          }
        ++randomWalkRegionIterator;
        }

      ANNVectorsCreated += currentBufferCounter;
      // std::cout<<" vec * curr :: "<<vector_size*currentBufferCounter<<"\n";
      ANNVectorStream.write( (const char *)( write_buffer ),
                             sizeof( float ) * vector_size * currentBufferCounter );
      delete[] write_buffer;
      }
      { /*Make an Error Image of improbable values*/
        // TODO Check where these images goes, and might have to change
        // writing location to the proper one. - KEY
      std::string ImprobableLocationsImageName = KnownOutputMaskName
        + "ImprobableImage.nii.gz";
      itkUtil::WriteImage<ProbabilityMapImageType>(ImprobableLocations,
                                                   ImprobableLocationsImageName);
      std::string DefProbImageName = KnownOutputMaskName
        + "DeformedProbMap.nii.gz";
      itkUtil::WriteImage<ProbabilityMapImageType>(DeformedProbabilityMap[currStructureIndex],
                                                   DefProbImageName);
      std::string DefRhoMapName = KnownOutputMaskName + "DeformedRhoMap.nii.gz";
      itkUtil::WriteImage<RealImageType>(DeformedRhoMap, DefRhoMapName);
      std::string DefPhiMapName = KnownOutputMaskName + "DeformedPhiMap.nii.gz";
      itkUtil::WriteImage<RealImageType>(DeformedPhiMap, DefPhiMapName);
      std::string DefThetaMapName = KnownOutputMaskName
        + "DeformedThetaMap.nii.gz";
      itkUtil::WriteImage<RealImageType>(DeformedThetaMap, DefThetaMapName);
      }
    }
  if( doSVM == true )
    {
    int SVMVectorIncrement = 1;
    if( ( NonZeroProbMapPoints > SVMTrainingSize ) && ( doSVM ) )
      {
      SVMVectorIncrement = ( NonZeroProbMapPoints ) / SVMTrainingSize;
      SVMVectorIncrement = ( SVMVectorIncrement > 0 ) ? SVMVectorIncrement : 1;
      std::cout << "(Non-zero prob-map/SVMTrainingSize): ("
                << NonZeroProbMapPoints << "/" << SVMTrainingSize << ")="
                << SVMVectorIncrement << "  " << min << "  " << max
                << std::endl;
      }
    int reducer = SVMVectorIncrement;
    bbri = bbri.Begin();
    while( !bbri.IsAtEnd() )
      {
      const ProbabilityMapImageType::IndexType CurrentIndex = bbri.GetIndex();
      ProbabilityMapImageType::PixelType       current_value =
        static_cast<ProbabilityMapImageType::PixelType>(
          DeformedProbabilityMap[currStructureIndex]
          ->GetPixel(
            CurrentIndex) );
      reducer--;
      if( ( current_value <
            ( HUNDRED_PERCENT_VALUE - PERCENT_MIN_MAX_TOLERANCE ) )
          &&  ( current_value > ( 0.0F + PERCENT_MIN_MAX_TOLERANCE ) )
          && ( reducer == 0 )
          )
        {
        // Write to the output stream
        reducer = SVMVectorIncrement;
        for( std::vector<float>::const_iterator ofi = outputvector.begin();
             ofi != outputvector.end();
             ++ofi )
          {
          if( *ofi == 0 )
            {
            SVMVectorStream << "-1 ";
            }
          else
            {
            SVMVectorStream << "+1 ";
            }
          }
        SVMVectorStream << "\t";
        for( std::vector<float>::const_iterator ifi = inputvector.begin();
             ifi != inputvector.end();
             ++ifi )
          {
          SVMVectorStream << *ifi << ' ';
          }
        SVMVectorStream << std::endl;
        SVMVectorsCreated++;
        }
      ++bbri;
      }
    }
  return true;
}

// / Create Training Vectors
/**
 * Creates Training Vectors For All Images
 * \ingroup Step2
 * \param ANNXMLObject AS Description
 */
int CreateVectors(ProcessDescription & ANNXMLObject,
                  bool histogramEqualization,
                  bool doTest,
                  int verbose)
{
  try
    {
    RegistrationParams *regParams =
      ANNXMLObject.Get<RegistrationParams>("RegistrationParams");
    if( regParams == 0 )
      {
      std::cout << "Registration information (parameter map filename, atlas image filename) "
                << "not given, so training vectors cannot be computed."
                << std::endl;
      exit(-1);
      }
    // Look for NeuralNetParams
    NeuralParams *model = ANNXMLObject.Get<NeuralParams>("NeuralNetParams");
    if( model == 0 )
      {
      std::cout << "NeurlaNetParams (gaussian size, iris size, "
                << "mask smoothing value, training vector filename) not given,"
                << " so training vectors cannot be created."
                << std::endl;
      exit(-1);
      }

    std::cout << "Using Training Filename file: "
              <<   model->GetAttribute<StringValue>("TrainingVectorFilename")
              << "\nUsing Mask Smoothing Value: "
              <<   model->GetAttribute<FloatValue>("MaskSmoothingValue")
              << "\nUsing GradientProfile Size: " << model->GetAttribute<IntValue>(
      "GradientProfileSize")
              << std::endl;

    bool doANN = 1;
    // Look for ANNParams
    ANNParams *annParams = ANNXMLObject.Get<ANNParams>("ANNParams");
    if( annParams == 0 )
      {
      doANN = 0;
      }

    bool doSVM = 1;
    int  SVMVectorSize = 0;
    // Look for SVMParams
    SVMParams *svModel = ANNXMLObject.Get<SVMParams>("SVMParams");
    //  if(svi == ANNXMLObject.end())
    if( svModel == 0 )
      {
      doSVM = 0;
      }
    else
      {
      SVMVectorSize = svModel->GetAttribute<IntValue>("VectorSize");
      }

    if( ( !doANN ) && ( !doSVM ) )
      {
      std::cout << "No models (ANN, SVM) chosen" << std::cout;
      return -1;
      }

    DataSet::TypeVector ImageTypeList =
      ANNXMLObject.GetAtlasDataSet()->ImageTypes();
    const int NumberOfImageTypes = ImageTypeList.size();
    if( !NumberOfImageTypes )
      {
      std::cout
        << "No images types found. Cannot compute neural net output." << std::endl;
      exit(-1);
      }

    std::cout << "------ "
              << model->GetAttribute<IntValue>("GradientProfileSize")
              << std::endl;

    std::string   ANNVectorFilename = "InvalidANNVectors.txt";
    std::ofstream ANNVectorStream;
    unsigned int  ANNTrainingVectorsCreated = 0;
    if( doANN )
      {
      if( !doTest )
        {
        ANNVectorFilename = model->GetAttribute<StringValue>("TrainingVectorFilename");
        }
      else if( doTest )
        {
        std::cout << " ******************************************************\n"
                  << " This is testing vector file generation \n"
                  << " Using Apply data set \n"
                  << " ******************************************************"
                  << std::endl;

        ANNVectorFilename = model->GetAttribute<StringValue>("TestVectorFilename");
        }
      ANNVectorFilename += "UnshuffledANN";

      std::string destination_dir =   itksys::SystemTools::GetFilenamePath(ANNVectorFilename);
      itksys::SystemTools::MakeDirectory( destination_dir.c_str() );
      ANNVectorStream.open(ANNVectorFilename.c_str(), std::ios::out | std::ios::binary);
      if( !ANNVectorStream.good() )
        {
        std::cout << "Error: Could not open ANN vector file: "
                  << ANNVectorFilename << std::endl;
        return -1;
        }
      }

    const std::string landmarkType("");
    const std::string SubjectLandmark("");
    const std::string AtlasLandmark("");

    const unsigned int numThreads = 1;

    // itk::MultiThreader::GetGlobalDefaultNumberOfThreads();
    GenerateRegistrations(ANNXMLObject, true, false, numThreads);

    const std::string regID( regParams->GetAttribute<StringValue>("ID") );

    const std::string imageTypeToUse
    (
      regParams->GetAttribute<StringValue>("ImageTypeToUse")
    );
    DataSet *     atlasDataSet = ANNXMLObject.GetAtlasDataSet();
    std::ofstream SVMVectorStream;
    unsigned int  SVMTrainingVectorsCreated = 0;
    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // Process 1. Getting  Output Vector Size == Number of ROIs
    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ProbabilityMapList *probabilityMaps =
      ANNXMLObject.Get<ProbabilityMapList>("ProbabilityMapList");
    int NumberOfProbabilityMaps = probabilityMaps->size();
    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // Process 2. Setting Input Vector Size Requirement
    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    const unsigned int inputVectorSize =
      InputVectorSizeRequirement(NumberOfProbabilityMaps,
                                 NumberOfImageTypes,
                                 model->GetAttribute<IntValue>("GradientProfileSize")
                                 );
    if( inputVectorSize == 0 )
      {
      std::cout << "Input vector size must not be zero!" << std::endl;
      exit(-1);
      }

    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // Process 3.
    //  -- Creating Hash Table for Each Prob Map
    //  -- Deform every Prob. Maps onto the Subject Space.
    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // In case of Testing Vector Creation, it has to use
    // 1. a TestFilename to create vector
    // 2. apply data set to create vector
    // Scan Each Target and Apply Trained Model to the Target
    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    ProcessDescription::TrainDataSetListType inputVectorDataSets;
    if( !doTest )
      {
      inputVectorDataSets = ProcessDescription::TrainDataSetListType(
          ANNXMLObject.GetTrainDataSets() );
      }
    else
      {
      std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n"
                << " This input vector creation is for test vector set. \n"
                << " The vector will be created using apply data set \n"
                << " instead of training data set\n"
                << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
                << std::endl;
      inputVectorDataSets = ProcessDescription::TrainDataSetListType(
          ANNXMLObject.GetApplyDataSets() );
      }
    /*ProcessDescription::TrainDataSetListType inputVectorDataSets(
                                        ANNXMLObject.GetTrainDataSets() );*/

    int DataSetCount = 0;
    for( ProcessDescription::TrainDataSetListType::iterator dsIt =
           inputVectorDataSets.begin();
         dsIt != inputVectorDataSets.end(); ++dsIt, ++DataSetCount )
      {
      const std::string       ImageID( ( *dsIt )->GetAttribute<StringValue>("Name") );
      const RegistrationType *transform = ( *dsIt )->GetRegistrationWithID(regID);
      const std::string       AtlasToSubjRegistrationFilename
      (
        transform->GetAttribute<StringValue>("AtlasToSubjRegistrationFilename")
      );

      std::cout << "Processing Dataset "
                << ( *dsIt )->GetAttribute<StringValue>("Name") << std::endl;
      // Reading Input Images to the Map
      std::map<std::string, ProbabilityMapImageType::Pointer> MapOfImages;
      for( DataSet::TypeVector::const_iterator ImageList = ImageTypeList.begin();
           ImageList != ImageTypeList.end();
           ++ImageList )
        {
        const std::string
          curFilename( ( *dsIt )->GetImageFilenameByType(*ImageList) );
        // HACK:  radius is hardcoded to size 2, it should be taken from command
        // line arguments.
        ProbabilityMapImageType::SizeType radius;
        radius[0] = 0; radius[1] = 0; radius[2] = 0;
        MapOfImages[*ImageList] =
          ReadMedianFilteredImage<ProbabilityMapImageType>(curFilename, radius);
        }
      // Regina::: Checking origin should go here
      for( DataSet::TypeVector::const_iterator ImageList = ImageTypeList.begin();
           ImageList != ImageTypeList.end();
           ++ImageList )
        {
        std::cout << "Check Image Origin" << std::endl;
        const std::string
          curFilename1( ( *dsIt )->GetImageFilenameByType(*ImageList) );
        if( ImageList != ImageTypeList.begin() )
          {
          const std::string
            curFilename2( ( *dsIt )->GetImageFilenameByType( *( ImageList - 1 ) ) );

          if( MapOfImages[*ImageList]->GetOrigin()
              != MapOfImages[*( ImageList - 1 )]->GetOrigin() )
            {
            std::cerr << " * ERROR : Input images are not in the same space! " << std::endl
                      << " * ERROR : This will cause bad registration!!!!    " << std::endl
                      << " * ERROR : File 1 :: " << MapOfImages[*ImageList]->GetOrigin()
                      << curFilename1
                      << " * ERROR : File 2 :: " << MapOfImages[*( ImageList + 1 )]->GetOrigin()
                      << curFilename2 << std::endl;
            }
          }
        }

      // Use first image in list to get spacing and size
      ProbabilityMapImageType::Pointer           ReferenceImage = MapOfImages[*( ImageTypeList.begin() )];
      const ProbabilityMapImageType::SpacingType ImageSpacing =
        ReferenceImage->GetSpacing();
      const ProbabilityMapImageType::PointType ImageOrigin =
        ReferenceImage->GetOrigin();
      const ProbabilityMapImageType::SizeType ImageSize =
        ReferenceImage->GetLargestPossibleRegion().GetSize();
      ProbabilityMapImageType::Pointer
          DeformedProbabilityMap[NumberOfProbabilityMaps];
      int tempIndex = 0;
      // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      // Initilization for For Loop
      // ProbabilityMapIndexType probabilityMapIndex;
      for(  ProbabilityMapList::iterator pmi = probabilityMaps->begin();
            pmi != probabilityMaps->end();
            ++pmi,  tempIndex++ )
        {
        ProbabilityMap *probabilityMap =
          dynamic_cast<ProbabilityMap *>( pmi->second );

        // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        //  -- Deform every Prob. Maps onto the Subject Space.

        std::cout << "Deforming "
                  << probabilityMap->GetAttribute<StringValue>("Filename")
                  << " into coordinate space of "
                  << ImageID
                  << " with " << AtlasToSubjRegistrationFilename
                  << std::endl;
        DeformedProbabilityMap[tempIndex] =
          ImageWarper<RealImageType>(AtlasToSubjRegistrationFilename,
                                     probabilityMap->GetAttribute<StringValue>("Filename"),
                                     ReferenceImage);
        if( verbose )
          {
          std::cout << "This is end of Deforming ..... " << std::endl;
          }
        } // END OF for (  ProbabilityMapList::iterator pmi =
          // probabilityMaps->begin();

      // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      //  -- Create Vectors For each Structure

      // * Now we can generate Input vector with all deformed prob. maps *
      tempIndex = 0;
      for(  ProbabilityMapList::iterator pmi = probabilityMaps->begin();
            pmi != probabilityMaps->end();
            ++pmi,  tempIndex++ )
        {
        ProbabilityMap *probabilityMap =
          dynamic_cast<ProbabilityMap *>( pmi->second );

        // start of do only for true of probability map
        if( probabilityMap->GetAttribute<StringValue>("GenerateVector") == "true" )
          {
          std::cout << "Generate Vectors for =======================================" << std::endl
                    << ImageID  <<  std::endl
                    << probabilityMap->GetAttribute<StringValue>("StructureID") << std::endl
                    << "============================================================" << std::endl;

          const std::string atlasImage( atlasDataSet->GetImageFilenameByType(
                                          imageTypeToUse.c_str() ) );
          bool dataset = GenerateInputOutputVectors
            (

              // ANNXMLObject,

              ( *dsIt ),
              probabilityMap->GetAttribute<StringValue>("Filename"),
              probabilityMap->GetAttribute<StringValue>("rho"),

              // rho
              probabilityMap->GetAttribute<StringValue>("phi"),

              // phi
              probabilityMap->GetAttribute<StringValue>("theta"),

              // theta
              probabilityMap->GetAttribute<StringValue>("StructureID"),
              DeformedProbabilityMap,
              atlasDataSet, // KEY_HE
                            // atlasImage, AtlasFilename imageTypeToUse,
                            // AtlasLandmark,
                            // probabilityMap->GetAttribute<FloatValue>("Gaussian"),
              ImageTypeList,
              AtlasToSubjRegistrationFilename,
              inputVectorSize,
              model->GetAttribute<IntValue>("GradientProfileSize"),
              model->GetAttribute<FloatValue>("MaskSmoothingValue"),
              doANN,
              // ANNVectorSize,
              ANNTrainingVectorsCreated,
              ANNVectorStream,
              false,
              SVMVectorSize,
              SVMTrainingVectorsCreated,
              SVMVectorStream,
              tempIndex,
              NumberOfProbabilityMaps,
              histogramEqualization,
              verbose
            );
          if( dataset == false )
            {
            std::cout << "  Fail to Generate Vector At this state.... " << std::endl;
            }
          }
        else
          {
          std::cout << "Skipping Generate Vectors for ====================================== " << std::endl
                    << ImageID  <<  std::endl
                    << probabilityMap->GetAttribute<StringValue>("StructureID") << std::endl
                    << "============================================================" << std::endl;
          }
        // end of do only for true of probability map
        } // END OF Second for (  ProbabilityMapList::iterator pmi =
          // probabilityMaps->begin();
      }   // END OF for ( ProcessDescription::TrainDataSetListType::iterator
          // dsIt

    /*
          if ( doSVM )
            {
            const std::string SVMHeaderVectorStreamHeaderFilename
              = SVMVectorFilename + ".hdr";
            std::ofstream SVMHeaderVectorStream;
            SVMHeaderVectorStream.open(
              SVMHeaderVectorStreamHeaderFilename.c_str(), std::ios::out
              | std::ios::binary );
            if ( !SVMHeaderVectorStream.good() )
              {
              std::cout << "Error: Could not open SVM vector file:"
                        << SVMHeaderVectorStreamHeaderFilename << std::endl;
              return -1;
              }

            SVMHeaderVectorStream
                                   << "IVS "
                                   << inputVectorSize
                                   << std::endl
                                   << "TVC "
                                   << static_cast<unsigned int>(
              SVMTrainingVectorsCreated )
                                   << std::endl;
            SVMHeaderVectorStream.close();
            }
         //END OF EACH PROBABILITY MAP SCANNING*/
    if( doANN )
      {
      ANNVectorStream.close();
      const std::string ANNHeaderVectorStreamHeaderFilename =
        ANNVectorFilename + ".hdr";
      std::ofstream ANNHeaderVectorStream;
      ANNHeaderVectorStream.open(
        ANNHeaderVectorStreamHeaderFilename.c_str(), std::ios::out
        | std::ios::binary);
      if( !ANNHeaderVectorStream.good() )
        {
        std::cout << "Error: Could not open ANN vector file: "
                  << ANNHeaderVectorStreamHeaderFilename << std::endl;
        return -1;
        }

      ANNHeaderVectorStream <<  "IVS " <<   inputVectorSize <<     std::endl
                            <<  "OVS " <<   NumberOfProbabilityMaps   <<    std::endl
                            <<   "TVC "
                            <<   static_cast<unsigned int>(
        ANNTrainingVectorsCreated )
                            <<     std::endl;
      ANNHeaderVectorStream.close();
      }
    }
  catch( ProcessObjectException & ex )
    {
    std::cerr << "Error1:: " << ex.Error() << std::endl;
    exit(-1);
    }
  return 0;
}

// for each probability map
//   for each image
//     create atlas-to-subject transform deformation field if needed
//     deform probability map into subject space
//     find mask for current image + structure id

int CreateVectors(const std::string & XMLFile,
                  bool histogramEqualization,
                  bool doTest,
                  int verbose)
{
  int rval = 0;

  try
    {
    ProcessDescription ANNXMLObject;

    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // + Read XML file and check its validity.
    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    try
      {
      ReadXML(XMLFile.c_str(), ANNXMLObject);
      }
    catch( ProcessObjectException & ex )
      {
      std::cerr << ex.Error() << std::endl;
      exit(-1);
      }
    if( ANNXMLObject.Verify() != true )
      {
      std::cerr << "XML file " << " is invalid." << std::endl;
      std::cerr << "FULL DUMP ===============================" << std::endl;
      std::cerr << "FULL DUMP ===============================" << std::endl;
      std::cerr << "FULL DUMP ===============================" << std::endl;
      ANNXMLObject.PrintSelf(std::cerr, 0);
      std::cerr << "FULL DUMP ===============================" << std::endl;
      std::cerr << "FULL DUMP ===============================" << std::endl;
      std::cerr << "FULL DUMP ===============================" << std::endl;
      return -1;
      }
    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // + Create Vectors
    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    try
      {
      rval = CreateVectors(ANNXMLObject,
                           histogramEqualization,
                           doTest,
                           verbose);
      }
    catch( ProcessObjectException & ex )
      {
      std::cerr << ex.Error() << std::endl;
      exit(-1);
      }
    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // + Shuffle Vectors
    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    try
      {
      NeuralParams *model = ANNXMLObject.Get<NeuralParams>("NeuralNetParams");

      // TFN = Trianing Filename
      std::string base = model->GetAttribute<StringValue>("TrainingVectorFilename");

      std::string TFN = base + "ANN";
      std::string unshuffledTFN = base + "UnshuffledANN";

      ShuffleVectors * my_ShuffleVector
        = new ShuffleVectors(  unshuffledTFN, TFN );
      my_ShuffleVector->ReadHeader();
      my_ShuffleVector->Shuffling();
      my_ShuffleVector->WriteHeader();
      }
    catch( ProcessObjectException & ex )
      {
      std::cerr << ex.Error() << std::endl;
      exit(-1);
      }
    rval = 0;
    }
  catch( ProcessObjectException & ex )
    {
    std::cerr << ex.Error() << std::endl;
    exit(-1);
    }
  catch( ... )
    {
    std::cerr << " Unknown execption while creating vectors "
              << std::endl;
    }
  return rval;
}
