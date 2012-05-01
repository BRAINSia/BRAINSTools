/**
 * \defgroup Step3 Step 3 - Train Neural Net
 * \ingroup Main
 * In the third step of the software, the neural network is trained.
 */
#include "Utilities.h"
#include "itksys/SystemTools.hxx"
#include <NetConfiguration.h>
#include "TrainingVectorConfigurationType.h"
#include "ANNParams.h"
// #include "ApplyModel.h"
#include "NetConfigurationParser.h"
#include <vnl/vnl_random.h>
#include <sys/time.h>

#if defined( NUMERIC_TRAP_TEST )
#include <xmmintrin.h>
#endif
#define EPS_MSE 1e-14F

unsigned long * ShuffledOrder(unsigned long shuffleSize)
{
  vnl_random     randgen;
  unsigned long *rval = new unsigned long[shuffleSize];

  if( rval == 0 )
    {
    itkGenericExceptionMacro(<< "Can't allocate shuffled ordering");
    }
  for( unsigned long i = 0; i < shuffleSize; i++ )
    {
    rval[i] = static_cast<std::ios::off_type>( i );
    }
  // do the shuffle
  for( unsigned long i = shuffleSize - 1; i > 0; i-- )
    {
    unsigned long j( randgen.lrand32() % ( i + 1 ) ); // rand() % (i+1);
    unsigned long tmp(rval[i]);
    rval[i] = rval[j];
    rval[j] = tmp;
    }
  return rval;
}

static unsigned long *vector_order = 0;
static unsigned long  vector_index = 0;
static unsigned       start_iteration = 0;

/** \author hjohnson
 * It is often impossible to read in an entire training vector complement at once.
 * This function will read treat the total set as having NumSubSets interleaved,
 * and will put out SubSetNumber.
 * For example if NumSubSets=4 and SubSetNumber=3, then elements 3,7,11,15...
 * would be read in.  If NumSubSets=4 and SubSetNumber=0, then elements
 * 4,8,12,16 would be used.
 */
neural_data_set_type *
ExtractTrainingSubsetFromFile(std::ifstream & vectorstr,  // const std::string
                                                          // ANNVectorFilename,
                              const int InputVectorSize,
                              const int OutputVectorSize,
                              const int NumSubSets,
                              const int SubSetNumber,
                              const int NumberTrainingVectorsFromFile,
                              int verbose = 0)
{
  int                   SliceSize    = NumberTrainingVectorsFromFile / NumSubSets;
  neural_data_set_type *TrainSetPtr =    new neural_data_set_type;
  neural_scalar_type *  tempBuf = new neural_scalar_type[SliceSize * InputVectorSize];
  neural_scalar_type *  tempOutputBuf = new neural_scalar_type[SliceSize * OutputVectorSize];

// ******************************************************************************//
  if( SubSetNumber >  NumSubSets )
    {
    itkGenericExceptionMacro(<< " Current Subset number ," << SubSetNumber
                             << " can not exceed total number of subset, " << NumSubSets);
    }
  int                NumInputVectors;
  std::ios::off_type recordSize =
    ( InputVectorSize + OutputVectorSize + 1 ) * sizeof( neural_scalar_type);
  std::ios::off_type  seekval;
  neural_scalar_type* buf = new neural_scalar_type[InputVectorSize + OutputVectorSize + 1];
  // two conditions pertain -- either you've gotten your full count
  // of vectors, or you've run out of file. In which case, you're done.
  // the training
  for( NumInputVectors = 0;
       ( NumInputVectors < SliceSize )   &&   !vectorstr.eof();
       NumInputVectors++ )
    {
    seekval = static_cast<std::ios::off_type>( vector_order[vector_index] )
      * recordSize;
    // vector_index = (vector_index ) % NumberTrainingVectorsFromFile;
    vector_index = ( vector_index + 1 ) % NumberTrainingVectorsFromFile;
    vectorstr.seekg(seekval, std::ios::beg);
    vectorstr.read( (char *)buf, recordSize );
    for( int i = 0; i < OutputVectorSize; i++ )
      {
      tempOutputBuf[NumInputVectors * OutputVectorSize + i] = buf[i];
      }
    for( int i = 0; i < InputVectorSize; i++ )
      {
      tempBuf[NumInputVectors * InputVectorSize + i] = buf[i + OutputVectorSize];
      }
    if( buf[InputVectorSize + OutputVectorSize] != AUTOSEG_VEC_SENTINEL )
      {
      itkGenericExceptionMacro(<< "Vector size mismatch in Vector file "
                               << "while reading vector number:"
                               << NumInputVectors << std::endl
                               << "buffer contents:: "
                               << buf[InputVectorSize + OutputVectorSize]);
      }
    }

  delete[] buf;
  std::cout << "\nUsing "
            << NumInputVectors
            << " for this iteration of training."
            << std::endl << std::flush;
  TrainSetPtr->inputVector = cvCreateMat( SliceSize, InputVectorSize, CV_32FC1);
  TrainSetPtr->outputVector = cvCreateMat( SliceSize, OutputVectorSize, CV_32FC1);
  cvInitMatHeader( TrainSetPtr->inputVector,
                   SliceSize,
                   InputVectorSize,
                   CV_32FC1,
                   tempBuf);
  cvInitMatHeader( TrainSetPtr->outputVector,
                   SliceSize,
                   OutputVectorSize,
                   CV_32FC1,
                   tempOutputBuf);

  if( verbose > 7 )
    {
    for( int i = 0; i < SliceSize; i++ )
      {
      std::cout << "\n I  ";
      for( int j = 0; j < InputVectorSize; j++ )
        {
        std::cout << " : " << CV_MAT_ELEM(  *(TrainSetPtr->inputVector), neural_scalar_type, i, j);
        }
      std::cout << " :: O";
      for( int j = 0; j < OutputVectorSize; j++ )
        {
        std::cout << " : " << CV_MAT_ELEM( *(TrainSetPtr->outputVector), neural_scalar_type, i, j);
        }
      }
    }
  return TrainSetPtr;
}

void ANNTrain( // NetConfiguration & prob,
  ANNParams *annParams,
  //  const std::string &StructureID,
  const std::string & VectorFilename,
  const std::string & ModelFilename,
  const std::string & TestVectorFilename,
  int verbose)
{
  std::string ANNVectorFilename = VectorFilename;

  ANNVectorFilename += "ANN";
  std::string ANNTestVectorFilename = TestVectorFilename;
  ANNTestVectorFilename += "ANN";
  std::string ANNModelFilename = ModelFilename;
  ANNModelFilename += "ANN";
  std::string ANNErrorFilename = ModelFilename;
  ANNErrorFilename += "ANNError";

  int TotalTrainingIterations =    annParams->GetAttribute<IntValue>("Iterations");
  // New Parameters
  int         MaximumVectorsPerEpoch =    annParams->GetAttribute<IntValue>("MaximumVectorsPerEpoch");
  const int   EpochIterations  =    annParams->GetAttribute<IntValue>("EpochIterations");
  const int   ErrorInterval    =    annParams->GetAttribute<IntValue>("ErrorInterval");
  const float DesiredError   =    annParams->GetAttribute<FloatValue>("DesiredError");
  int         HiddenVectorSize =    annParams->GetAttribute<IntValue>("NumberOfHiddenNodes");
  /*
   * openCV parameter
   *   a = alpha = slope = sigmoid_slope
   *   b = beta = min/max = sigmoid_range
   *   --> f(x) = b * ( 1 - exp( -ax ) / (1 + exp( -ax ) )
   */
  const float sigmoid_slope =    annParams->GetAttribute<FloatValue>("ActivationSlope");
  const float sigmoid_range =    annParams->GetAttribute<FloatValue>("ActivationMinMax");

  int               NumberTrainingVectorsFromFile = 0;
  int               InputVectorSize = 0;
  int               OutputVectorSize = 0;
  const std::string ANNHeaderVectorFilename = ANNVectorFilename + ".hdr";
  std::ifstream     filestr;
  filestr.open(ANNHeaderVectorFilename.c_str(), std::ios::in | std::ios::binary);
  if( !filestr.is_open() )
    {
    itkGenericExceptionMacro(<< "Error:  ANNTrain: Could not open ANN vector file name of \""
                             << ANNHeaderVectorFilename.c_str() << "\"");
    }
  for( int tags = 0; tags < 3; tags++ )  // Read header file for input/output
                                         // vectors and total vectors.
    {
    std::string temp;
    char        currentline[MAX_LINE_SIZE];
    filestr.getline(currentline, MAX_LINE_SIZE - 1);
    std::istringstream iss(currentline, std::istringstream::in);
    iss >> temp;
    if( temp == "IVS" )
      {
      iss >> InputVectorSize;
      }
    if( temp == "OVS" )
      {
      iss >> OutputVectorSize;
      }
    else if( temp == "TVC" )
      {
      iss >> NumberTrainingVectorsFromFile;
      }
    }
  filestr.close();
  if( ( InputVectorSize == 0 ) || ( NumberTrainingVectorsFromFile == 0 )
      || ( OutputVectorSize == 0 ) )
    {
    std::cout
      << "Input/output vector size and vector count must be non-zero."
      << std::endl;
    return;
    }
  // set up shuffled ordering to read out of training vector file
  vector_order = ShuffledOrder(NumberTrainingVectorsFromFile);
  if( HiddenVectorSize == 0 )
    {
    HiddenVectorSize =  ( InputVectorSize ) * 2;
    }
  std::cout << "\nNeural Net Configurations:"
            << "\nNet Weights Filename: "  << ANNModelFilename
            << "\nTraining Filename: " << ANNVectorFilename
            << "\nError Filename: " << ANNErrorFilename
            << "\nTraining Iterations: " << TotalTrainingIterations
            << "\nError Interval: " << ErrorInterval
            << "\nNumber Training Vectors Per Error Interval: "
            << "\nTotal Number of Vectors: "
            << NumberTrainingVectorsFromFile << std::endl;
  std::cout << "\nEstimated Memory Needed for Vectors: (MB) "
            << ( ( static_cast<float>( MaximumVectorsPerEpoch ) )
       * sizeof( neural_scalar_type)
       * ( InputVectorSize + 1.0F ) ) / ( 1024.0F * 1024.0F )
            << std::endl;

  // Provide a random number seed value so that multiple identical runs produce
  // the same result.
  srand(123);

  std::string Weights_Filename(ANNModelFilename);
  if( start_iteration != 0 )
    {
    char num[32];
    sprintf(num, "%09u", start_iteration);
    Weights_Filename += num;
    if( !itksys::SystemTools::FileExists( Weights_Filename.c_str() ) )
      {
      std::cerr << "Weights File " << Weights_Filename
                << " does not exist." << std::endl;
      throw;
      }
    std::cout << "Using initilization File name of .. "
              << Weights_Filename << std::endl;
    }
  std::cout << "\nTraining Neural Net..." << std::endl;
  std::ofstream ErrorStream;
    {
    std::string destination_dir = itksys::SystemTools::GetFilenamePath(
        ANNErrorFilename);
    itksys::SystemTools::MakeDirectory( destination_dir.c_str() );
    }
  ErrorStream.open(ANNErrorFilename.c_str(), std::ios::out | std::ios::binary);
  if( !ErrorStream.good() )
    {
    std::cout << "Could Not Open Error File: "  << ANNErrorFilename
              << std::endl;
    throw;
    }

  struct timeval t_start, t_current;
  int            t_diff = 0;
  gettimeofday(&t_start, NULL);
  // HACK  This should be a command line variable, or related to the number of
  // epochs
  int NumSubSets = 0;
  if( MaximumVectorsPerEpoch == 0 )
    {
    std::cerr << "MaximumVectorsPerEpoch cannot be zero, "
              << "please check your xml file.\n";
    }
  else
    {
    NumSubSets =  NumberTrainingVectorsFromFile / MaximumVectorsPerEpoch + 1;
    }
  if( NumSubSets == 0 )
    {
    NumSubSets  = 1;
    }
  ErrorStream << "Mean Square Error for Neural Net:" << std::endl;

  int CurrentTimesTrained = start_iteration + ErrorInterval;

  // open training file once per feature here, instead of doing
  // it every time you read in a subset...
  std::ifstream vectorstr;
  vectorstr.open(ANNVectorFilename.c_str(), std::ios::in);

  int                   lastSubSetNumber = -1;
  neural_data_set_type *TrainSetPtr = 0;
  // const std::string     ANNHeaderTestVectorFilename = ANNTestVectorFilename + ".hdr";
  std::ifstream testfilestr;
  // TODO: DELETE neural_data_set_type *TestSetPtr = 0;
  // TODO: DELETE int                   test_InputVectorSize = 0;
  // TODO: DELETE int                   test_OutputVectorSize = 0;
  // TODO: DELETE int                   test_NumberTrainingVectorsFromFile = 0;
  std::ifstream testVectorStr;
  testfilestr.close();
  // To trace minimum traininig point
  double TrainSet_MinimumMSE = 100.0; int TrainSet_MinimumMSEPoint = 0;

  // To trace minimum traininig point
  // TODO: DELETE double TestSet_MinimumMSE = 100.0;
  // TODO: DELETE int TestSet_MinimumMSEPoint = 0;
  // OPENCV Training Model Creation
  int                layer[] = { InputVectorSize, HiddenVectorSize, OutputVectorSize };
  neural_vector_type layerStructure = cvCreateMat( 1, 3, CV_32SC1 );
  cvInitMatHeader( layerStructure, 1, 3, CV_32SC1, layer );

  neural_net_type * training_net = new neural_net_type();
  training_net->create( layerStructure,
                        CvANN_MLP::SIGMOID_SYM,
                        sigmoid_slope,
                        sigmoid_range);
  // Initial Training
  std::cout << "* Initial Training ... ... "
            << std::endl;
  TrainSetPtr = ExtractTrainingSubsetFromFile(vectorstr,
                                              InputVectorSize,
                                              OutputVectorSize,
                                              NumSubSets,
                                              0,
                                              NumberTrainingVectorsFromFile,
                                              verbose
                                              );
  int first_iter = training_net->train( (TrainSetPtr->inputVector),
                                        (TrainSetPtr->outputVector),
                                        NULL, // SAMPLE WEIGHTS
                                        0,    // sampel index
                                        CvANN_MLP_TrainParams( cvTermCriteria( CV_TERMCRIT_ITER
                                                                               | CV_TERMCRIT_EPS,
                                                                               EpochIterations, DesiredError),
                                                               CvANN_MLP_TrainParams::RPROP, //
                                                                                             // 3,
                                                               0.1,                          //
                                                                                             // _param1
                                                               FLT_EPSILON                   //
                                                                                             // _param2
                                                               ),
                                        ( 0
                                          // + CvANN_MLP::NO_INPUT_SCALE
                                          // +  CvANN_MLP::NO_OUTPUT_SCALE
                                          // flags, WEIGHT UPDATE=NO
                                        )
                                        );
  std::cout << " First iteration took " << first_iter << std::endl;
  std::string initialTrainingFile = ANNModelFilename + "000000001";
  std::cout << "Saving files of "
            << initialTrainingFile.c_str()
            << std::endl;
  training_net->save( initialTrainingFile.c_str() );
  for( int SubSetNumber = 0;
       CurrentTimesTrained < TotalTrainingIterations;
       ( CurrentTimesTrained += ErrorInterval ),
       ( lastSubSetNumber = SubSetNumber ),
       ( SubSetNumber =
           ( ( SubSetNumber == NumSubSets - 1 ) ? 0 : SubSetNumber + 1 ) ) )
    {
    double EpochError = 1.0;
    char   tempid[10];
    int    TempTrained = CurrentTimesTrained; // + ErrorInterval;

    sprintf(tempid, "%09d", TempTrained + 1 );
    std::string temp_Weights = ANNModelFilename + tempid;
    std::cout << " * Number of Train:   " << TempTrained << std::endl;
    std::cout << " * Number of CurrentTimesTrained: " << CurrentTimesTrained << std::endl;
    std::cout << " * File name ::" << temp_Weights.c_str() << std::endl;
    if( !itksys::SystemTools::FileExists( temp_Weights.c_str() ) )
      {
      // only read in new subset if it's different than the last one
      if( SubSetNumber != lastSubSetNumber )
        {
        // get rid of old subset
        if( TrainSetPtr != 0 )
          {
          delete TrainSetPtr;
          }
        std::cout << "Reading training vector subset (" << SubSetNumber + 1 << " of "
                  << NumSubSets << ") for training." << std::endl;
        TrainSetPtr = ExtractTrainingSubsetFromFile(vectorstr,
                                                    InputVectorSize,
                                                    OutputVectorSize,
                                                    NumSubSets,
                                                    SubSetNumber + 1,
                                                    NumberTrainingVectorsFromFile,
                                                    false
                                                    );
        }

      std::cout << "Training." << std::endl;

      std::cout << "training opencv\n";
      if( verbose > 8 )
        {
        int SliceSize    = NumberTrainingVectorsFromFile / NumSubSets;
        for( int i = 0; i < SliceSize; i++ )
          {
          std::cout << "\n I  [" << InputVectorSize << "] ";
          for( int j = 0; j < InputVectorSize; j++ )
            {
            std::cout << ": "
                      << CV_MAT_ELEM(  *(TrainSetPtr->inputVector), neural_scalar_type, i, j);
            }
          std::cout << " :: O [" << OutputVectorSize << "] ";
          for( int j = 0; j < OutputVectorSize; j++ )
            {
            std::cout << ": "
                      << CV_MAT_ELEM(  *(TrainSetPtr->outputVector), neural_scalar_type, i, j);
            }
          }
        }
      int iteration = training_net->train( (TrainSetPtr->inputVector),
                                           (TrainSetPtr->outputVector),
                                           NULL,   // SAMPLE WEIGHTS
                                           0,      // sampel index
                                           CvANN_MLP_TrainParams( cvTermCriteria
                                                                    ( CV_TERMCRIT_ITER
                                                                    | CV_TERMCRIT_EPS,
                                                                    EpochIterations,
                                                                    DesiredError),
                                                                  CvANN_MLP_TrainParams::RPROP, //
                                                                                                // 3,
                                                                  0.1,                          //
                                                                                                // _param1
                                                                  FLT_EPSILON                   //
                                                                                                // _param2
                                                                  ),
                                           ( CvANN_MLP::UPDATE_WEIGHTS
                                             // + CvANN_MLP::NO_INPUT_SCALE
                                             // +  CvANN_MLP::NO_OUTPUT_SCALE'
                                           ) // flags
                                           );

      training_net->save( temp_Weights.c_str() );

      std::cout << "Number of Iteration :: "
                << iteration << std::endl;

      EpochError = training_net->get_MSE();
      ErrorStream << "Iterations: "    << TempTrained
                  << "\tError: "       << EpochError
                  << "\tTime: "        << t_diff
                  << " ms"             << std::endl     << std::flush;
      std::cout << "Iterations: "    << TempTrained
                << "\tError: "       << EpochError
                << "\tTime: "        << t_diff
                << " ms"             << " (in minutes "      << t_diff
        / ( 1000.0F * 60 ) << " min)"  << std::endl  << std::flush;

      if( TrainSet_MinimumMSE > EpochError )
        {
        TrainSet_MinimumMSE = EpochError;
        TrainSet_MinimumMSEPoint = TempTrained;
        }
      gettimeofday(&t_current, NULL);
      t_diff =
        (int)( 1000.0 * ( t_current.tv_sec - t_start.tv_sec )
               + ( t_current.tv_usec - t_start.tv_usec ) / 1000.0 );
      }
    else
      {
      std::cout << "trained model already exist. " << std::endl;
      }
    /* Do Test is not working with OpenCV
#ifndef USE_OPENCV
    if ( doTest )
    {
      std::cout << "create model from :: "
                << temp_Weights.c_str()
                << std::endl;
      training_net.create_from_file( temp_Weights.c_str() );
      std::cout << "call testing..." << std::endl;
      training_net.test_data(*TestSetPtr);
      std::cout << "Iterations: "    << TempTrained
                << "\tPerformance Error: " << training_net.get_MSE()
                << std::endl;
      if( TestSet_MinimumMSE > training_net.get_MSE() ){
        TestSet_MinimumMSE =training_net.get_MSE();
        TestSet_MinimumMSEPoint = TempTrained;
      }
      // delete TestSetPtr;
    }
#endif */
    // Save Trainning Model

    if( EpochError < DesiredError )
      {
      std::cout << "Desired Error Has been achieved. Stop Training. \n";
      break;
      }
    }
  // Printing out minimum training error
  // TODO: Should I add this to the Header file??!!!! KEY
  if( TrainSet_MinimumMSEPoint != 0 )  // if training has been proceeeded
    {
    std::cout <<  "###################################################\n"
              << " * Minimum MSE of Training : " << TrainSet_MinimumMSE
              << " at "                          << TrainSet_MinimumMSEPoint
              << " Iteration. "
              << "\n###################################################\n";
    }
  // exit from loop bypasses deleting training set
  delete TrainSetPtr;

  ErrorStream.close();
  remove( ANNModelFilename.c_str() );
  std::cout
    << "\tANN Model "
    //    << StructureID
    << " written to file."
    << std::endl;
}

static int Train(NetConfiguration & ANNConfiguration, int verbose)
{
  TrainingVectorConfigurationType *model = ANNConfiguration.Get<TrainingVectorConfigurationType>(
      "TrainingVectorConfiguration");

  ANNParams *annParams = ANNConfiguration.Get<ANNParams>("ANNParams");

  if( annParams != 0 )
    {
    annParams = ANNConfiguration.Get<ANNParams>("ANNParams");
    ANNTrain(   // ANNConfiguration,
      annParams,
      model->GetAttribute<StringValue>("TrainingVectorFilename"),
      model->GetAttribute<StringValue>("TrainingModelFilename"),
      model->GetAttribute<StringValue>("TestVectorFileName"),
      verbose
      );
    }
  return 0;
}

int TrainModel(const std::string & XMLFile, int verbose, int StartIteration)
{
#if defined( NUMERIC_TRAP_TEST )
  // turn on numeric exceptions...
  _mm_setcsr( _MM_MASK_MASK & ~
              ( _MM_MASK_OVERFLOW | _MM_MASK_INVALID | _MM_MASK_DIV_ZERO ) );
#endif
  NetConfiguration *     ANNConfiguration;
  NetConfigurationParser ANNConfigurationParser = NetConfigurationParser( XMLFile );
  ANNConfiguration = ANNConfigurationParser.GetNetConfiguration();

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
  // restart training
  start_iteration = StartIteration;
  int rval;
  try
    {
    if( verbose )
      {
      std::cerr << "This will print out in detail with process....\n";
      }
    std::cout << "Training Model(s)" << std::endl;
    rval = Train(*ANNConfiguration, verbose);
    std::cout << "Finished Training" << std::endl;
    }
  catch( BRAINSCutExceptionStringHandler & ex )
    {
    std::cerr << ex.Error() << std::endl;
    rval = -1;
    }
  catch( ... )
    {
    std::cerr << "Unidentified exception" << std::endl;
    rval = -1;
    }
  return rval;
}
