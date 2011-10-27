#include <NetConfiguration.h>
#include <NeuralParams.h>
#include <ANNParams.h>
#include <itksys/SystemTools.hxx>
#include <list>
#include <vector>
#include <iostream>
#include <itkIO.h>
#include "Utilities.h"
#include "XMLParser.h"
#include "NetConfigurationParser.h"
#include "itkImageRandomConstIteratorWithIndex.h"
#include "itkHistogramMatchingImageFilter.h" // KEY_HE
#include "ShuffleVectors.h"
#include "itkLabelStatisticsImageFilter.h"

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
 */

#if 0 // TODO: DELETE
/**
 * Create Training Vector
 */
static int CreateVectorsFromNetConfiguration(NetConfiguration & ANNConfiguration,
                                             bool /* TODO:  DELETE THIS LINE histogramEqualization */,
                                             int /* TODO: DELTE THIS LINE verbose */ )
{
  /** Check Necessary Inputs */
  /** 1. registration params */
  RegistrationConfigurationParser *regParams =
    ANNConfiguration.Get<RegistrationConfigurationParser>("RegistrationConfigurationParser");

  if( regParams == 0 )
    {
    std::cout << "Registration information (parameter map filename, atlas image filename) "
              << "not given, so training vectors cannot be computed."
              << std::endl;
    exit(-1);
    }
  /** 2. neural network params */
  NeuralParams *model = ANNConfiguration.Get<NeuralParams>("NeuralNetParams");
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
            << "\nUsing GradientProfile Size: "
            << model->GetAttribute<IntValue>("GradientProfileSize")
            << std::endl;

  /** 3. ANN specific parameter*/
  ANNParams *annParams = ANNConfiguration.Get<ANNParams>("ANNParams");
  if( annParams == 0 )
    {
    std::cout << "No models ANN chosen" << std::cout;
    return -1;
    }

  /** 4. Get Image Types Ex. T1,T2 and SG*/
  DataSet::StringVectorType ImageTypeList =
    ANNConfiguration.GetAtlasDataSet()->GetImageTypes();
  const int NumberOfImageTypes = ImageTypeList.size();
  if( !NumberOfImageTypes )
    {
    std::cout
      << "No images types found. Cannot compute neural net output." << std::endl;
    exit(-1);
    }

  /** number of thread = 1 hard coded */
  const unsigned int numThreads = 1;

  /** generate registration from atlas to subject direction */
  GenerateRegistrations(ANNConfiguration, true, false,  numThreads);

  const std::string regID( regParams->GetAttribute<StringValue>("ID") );

  const std::string imageTypeToUse
  (
    regParams->GetAttribute<StringValue>("ImageTypeToUse")
  );

  /** atals */
  // TODO: Delete this DataSet * atlasDataSet = ANNConfiguration.GetAtlasDataSet();
  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // Process 1. Getting  Output Vector Size == Number of ROIs
  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ProbabilityMapList *probabilityMaps =
    ANNConfiguration.Get<ProbabilityMapList>("ProbabilityMapList");
  int NumberOfROIs = probabilityMaps->size();

  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // Process 2. Setting Input Vector Size Requirement
  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  const unsigned int inputVectorSize = InputVectorSizeRequirement(NumberOfROIs,
                                                                  NumberOfImageTypes,
                                                                  model->GetAttribute<IntValue>("GradientProfileSize") );
  if( inputVectorSize == 0 )
    {
    std::cout << "Input vector size must not be zero!" << std::endl;
    exit(-1);
    }

  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // Process
  //  -- Fixe the order of ROIs
  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ProbabilityMapList * ROIList = ANNConfiguration.Get<ProbabilityMapList>("ProbabilityMapList");

  std::map<int, std::string> MapOfROIOrder;
  int                        roiCounter = 0;
  for( ProbabilityMapList::iterator pmi = ROIList->begin();
       pmi != ROIList->end();
       ++pmi )
    {
    MapOfROIOrder.insert( pair<int, std::string>( roiCounter, pmi->first) );
    roiCounter++;
    }

  /** create file */
  std::string ANNVectorFilename = "InvalidANNVectors.txt";

  ANNVectorFilename = model->GetAttribute<StringValue>("TrainingVectorFilename");
  ANNVectorFilename += "UnshuffledANN";

  std::string destination_dir =   itksys::SystemTools::GetFilenamePath(ANNVectorFilename);
  itksys::SystemTools::MakeDirectory( destination_dir.c_str() );
  std::remove( ANNVectorFilename.c_str() );

  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // interate all the subjects
  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  NetConfiguration::TrainDataSetListType inputSubjectsList;

  inputSubjectsList = NetConfiguration::TrainDataSetListType(
      ANNConfiguration.GetTrainDataSets() );

  // Create vectors for each subject
  int SubjectCount = 0;
  int ANNTrainingVectorsCreated = 0;
  for( NetConfiguration::TrainDataSetListType::iterator subjectIt = inputSubjectsList.begin();
       subjectIt != inputSubjectsList.end();
       ++subjectIt, ++SubjectCount )
    {
    std::cout << __LINE__ << " :: " << __FILE__
              << "    :: " << *subjectIt << std::endl
              << "    :: " << ANNTrainingVectorsCreated << std::endl;
    ANNTrainingVectorsCreated += AddSubjectInputVector(
        (*subjectIt),
        ANNConfiguration,
        regID,
        inputVectorSize,
        NumberOfROIs,
        MapOfROIOrder,
        false);                      // last parameter: apply=false
    }

  /** ANNVectorStream */

  const std::string ANNHeaderVectorStreamHeaderFilename = ANNVectorFilename + ".hdr";
  std::ofstream     ANNHeaderVectorStream;
  ANNHeaderVectorStream.open( ANNHeaderVectorStreamHeaderFilename.c_str(),
                              std::ios::out | std::ios::binary);
  if( !ANNHeaderVectorStream.good() )
    {
    std::cout << "Error: Could not open ANN vector file: "
              << ANNHeaderVectorStreamHeaderFilename << std::endl;
    return -1;
    }

  ANNHeaderVectorStream <<  "IVS " <<   inputVectorSize <<     std::endl
                        <<  "OVS " <<   NumberOfROIs   <<    std::endl
                        <<   "TVC "
                        <<   static_cast<unsigned int>(  ANNTrainingVectorsCreated )
                        <<     std::endl;

  ANNHeaderVectorStream.close();
  return 0;
}

#endif

#if 0 // TODO:  DELETE
// for each probability map
//   for each image
//     create atlas-to-subject transform deformation field if needed
//     deform probability map into subject space
//     find mask for current image + structure id

static int CreateVectorsFromXMLFile(const std::string & XMLFile,
                                    bool /* TODO DELETE THIS histogramEqualization */,
                                    int /* verbose */ )
{
  int rval = 0;

  try
    {
    NetConfiguration *     ANNConfiguration;
    NetConfigurationParser ANNConfigurationParser = NetConfigurationParser( XMLFile );

    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // + Read XML file and check its validity.
    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    try
      {
      ANNConfiguration = ANNConfigurationParser.GetNetConfiguration();
      }
    catch( BRAINSCutExceptionStringHandler & ex )
      {
      std::cerr << ex.Error() << std::endl;
      exit(-1);
      }
    if( ANNConfiguration->Verify() != true )
      {
      std::cerr << "XML file " << " is invalid." << std::endl;
      std::cerr << "FULL DUMP ===============================" << std::endl;
      std::cerr << "FULL DUMP ===============================" << std::endl;
      std::cerr << "FULL DUMP ===============================" << std::endl;
      ANNConfiguration->PrintSelf(std::cerr, 0);
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
      rval = CreateVectorsFromNetConfiguration(*ANNConfiguration,
                                               true /* TODO DELETE THIS histogramEqualization */,
                                               true /* TODO DELETE THIS verbose */);
      }
    catch( BRAINSCutExceptionStringHandler & ex )
      {
      std::cerr << ex.Error() << std::endl;
      exit(-1);
      }
    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // + Shuffle Vectors
    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    try
      {
      NeuralParams *model = ANNConfiguration->Get<NeuralParams>("NeuralNetParams");

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
    catch( BRAINSCutExceptionStringHandler & ex )
      {
      std::cerr << ex.Error() << std::endl;
      exit(-1);
      }
    rval = 0;
    }
  catch( BRAINSCutExceptionStringHandler & ex )
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

#endif
