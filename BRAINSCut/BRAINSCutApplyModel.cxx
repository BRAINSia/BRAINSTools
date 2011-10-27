#include "BRAINSCutApplyModel.h"
#include "FeatureInputVector.h"
#include "ANNParams.h"
#include "Utilities.h"

BRAINSCutApplyModel
::BRAINSCutApplyModel( std::string netConfigurationFilename)
  : BRAINSCutPrimary( netConfigurationFilename )
  // ANNModelFilename(NULL)
{
  // TODO Take this apart to generate registration one by one!
  GenerateRegistrations(BRAINSCutNetConfiguration, true, true, 1);

  SetRegionsOfInterestFromNetConfiguration();
  SetRegistrationParametersFromNetConfiguration();
  SetAtlasDataSet();
  SetAtlasImage();
  SetRhoPhiThetaFromNetConfiguration();
  SetANNModelConfiguration();
  SetGradientSizeFromNetConfiguration();

  normalization = GetNormalizationFromNetConfiguration();

  openCVANN = new OpenCVMLPType();
}

/* iterate through subject */
void
BRAINSCutApplyModel
::Apply()
{
  typedef NetConfiguration::ApplyDataSetListType::iterator
    ApplySubjectIteratorType;
  for( ApplySubjectIteratorType subjectIt = applyDataSetList.begin();
       subjectIt != applyDataSetList.end();
       ++subjectIt )
    {
    ApplyOnSubject( *(*subjectIt) );
    }
}

/* apply on each subject */
void
BRAINSCutApplyModel
::ApplyOnSubject( DataSet& subject)
{
  std::map<std::string, WorkingImagePointer> deformedSpatialLocationImageList;

  GetDeformedSpatialLocationImages( deformedSpatialLocationImageList, subject );

  WorkingImageVectorType imagesOfInterest;
  GetImagesOfSubjectInOrder(imagesOfInterest, subject);

  std::map<std::string, WorkingImagePointer> deformedROIs;
  GetDeformedROIs(deformedROIs, subject);

  FeatureInputVector inputVectorGenerator;

  inputVectorGenerator.SetGradientSize( gradientSize );
  inputVectorGenerator.SetImagesOfInterestInOrder( imagesOfInterest );
  inputVectorGenerator.SetImagesOfSpatialLocation( deformedSpatialLocationImageList );
  inputVectorGenerator.SetCandidateROIs( deformedROIs);
  inputVectorGenerator.SetROIInOrder( roiIDsInOrder);
  inputVectorGenerator.SetInputVectorSize();
  inputVectorGenerator.SetNormalization( normalization );

  /* now iterate through the roi */

  unsigned int roiIDsOrderNumber = 0;
  for( DataSet::StringVectorType::iterator roiTyIt = roiIDsInOrder.begin();
       roiTyIt != roiIDsInOrder.end();
       ++roiTyIt ) // roiTyIt = Region of Interest Type Iterator
    {
    ProbabilityMapParser* roiDataSet =
      roiDataList->GetMatching<ProbabilityMapParser>( "StructureID", (*roiTyIt).c_str() );
    if( roiDataSet->GetAttribute<StringValue>("GenerateVector") == "true" )
      {
      InputVectorMapType  roiInputVector = inputVectorGenerator.GetFeatureInputOfROI( *roiTyIt );
      PredictValueMapType predictedOutputVector;
      PredictROI( roiInputVector, predictedOutputVector,
                  roiIDsOrderNumber, inputVectorGenerator.GetInputVectorSize() );

      std::string ANNContinuousOutputFilename = GetContinuousPredictionFilename( subject, (*roiTyIt) );

      WritePredictROIProbabilityBasedOnReferenceImage( predictedOutputVector,
                                                       imagesOfInterest.front(),
                                                       deformedROIs.find( *roiTyIt )->second,
                                                       ANNContinuousOutputFilename );
      /* post processing
       * may include hole-filling(closing), thresholding, and more adjustment
       */
      BinaryTypePointer mask = PostProcessingOfANNContinuousImage( ANNContinuousOutputFilename);

      std::string roiOutputFilename = GetOutputROIFilename( subject, *roiTyIt );
      itkUtil::WriteImage<BinaryImageType>( mask, roiOutputFilename );
      }

    // TODO:do clean up here
    // TODO:writing mask here ?
    roiIDsOrderNumber++;
    }
}

void
BRAINSCutApplyModel
::SetApplyDataSetFromNetConfiguration()
{
  try
    {
    applyDataSetList = BRAINSCutNetConfiguration.GetApplyDataSets();
    }
  catch( BRAINSCutExceptionStringHandler& e )
    {
    std::cout << e.Error() << std::endl;
    exit(EXIT_SUCCESS);
    }
}

BinaryTypePointer
BRAINSCutApplyModel
::PostProcessingOfANNContinuousImage( std::string continuousFIlename )
{
  WorkingImagePointer continuousImage = ReadImageByFilename( continuousFIlename );

  /* threshold first */
  BinaryTypePointer maskVolume;

  maskVolume = ThresholdImage( continuousImage );

  /* TODO hole filling here */

  return maskVolume;
}

BinaryTypePointer
BRAINSCutApplyModel
::ThresholdImage( WorkingImagePointer image )
{
  typedef itk::BinaryThresholdImageFilter<WorkingImageType, WorkingImageType> ThresholdFilterType;
  ThresholdFilterType::Pointer thresholder = ThresholdFilterType::New();

  thresholder->SetInput( image );
  thresholder->SetInsideValue( 1 );
  thresholder->SetOutsideValue( 0 );
  thresholder->SetLowerThreshold( annOutputThreshold  );
  thresholder->Update();

  BinaryTypePointer mask = itkUtil::TypeCast<WorkingImageType, BinaryImageType>( thresholder->GetOutput() );
  return mask;
}

void
BRAINSCutApplyModel
::SetANNOutputThresholdFromNetConfiguration()
{
  annOutputThreshold =
    BRAINSCutNetConfiguration.Get<ANNParams>("ANNParams")->GetAttribute<FloatValue>("MaskThresh");
  if( annOutputThreshold < 0.0F )
    {
    std::string msg = " ANNOutput Threshold cannot be less than zero. \n";
    throw BRAINSCutExceptionStringHandler( msg );
    }
}

inline void
BRAINSCutApplyModel
::PredictROI( InputVectorMapType&  roiInputFeatureVector,
              PredictValueMapType& resultOutputVector,
              unsigned int         roiNumber,
              unsigned int         inputVectorSize)
{
  ReadANNModelFile();
  for( InputVectorMapType::iterator it = roiInputFeatureVector.begin();
       it != roiInputFeatureVector.end();
       ++it )
    {
    /* convert input feature vector to array */
    scalarType* arrayInputFeatureVector = new scalarType[inputVectorSize];
    arrayInputFeatureVector = GetArrayFromVector( arrayInputFeatureVector,  it->second, inputVectorSize);

    /* get open cv type matrix from array for prediction */
    matrixType openCVInputFeature = cvCreateMat( 1, inputVectorSize, CV_32FC1);
    GetOpenCVMatrixFromArray( openCVInputFeature, arrayInputFeatureVector, inputVectorSize);

    /* predict */
    matrixType openCVOutput = cvCreateMat( 1, roiIDsInOrder.size(), CV_32FC1);
    openCVANN->predict( openCVInputFeature, openCVOutput );

    /* insert result to the result output vector */
    resultOutputVector.insert( pair<int, scalarType>(  ( it->first ),  CV_MAT_ELEM( *openCVOutput,
                                                                                    scalarType,
                                                                                    0,
                                                                                    roiNumber) ) );
    }
}

inline void
BRAINSCutApplyModel
::GetOpenCVMatrixFromArray( matrixType& matrix, scalarType array[], unsigned int inputVectorSize)
{
  cvInitMatHeader( matrix, 1, inputVectorSize, CV_32FC1, array );
}

void
BRAINSCutApplyModel
::SetANNModelFilenameFromNetConfiguration()
{
  try
    {
    ANNModelFilename = annModelConfiguration->GetAttribute<StringValue>("TrainingModelFilename");
    }
  catch( ... )
    {
    throw BRAINSCutExceptionStringHandler("Fail to get the ann model file name");
    }
  int  iteration = BRAINSCutNetConfiguration.Get<ANNParams>("ANNParams")->GetAttribute<IntValue>("Iterations");
  char temp[10];
  sprintf( temp, "%09d", iteration );
  ANNModelFilename += temp;
}

void
BRAINSCutApplyModel
::ReadANNModelFile()
{
  if( !itksys::SystemTools::FileExists( ANNModelFilename.c_str() ) )
    {
    std::string errorMsg = " File does not exist! :";
    errorMsg += ANNModelFilename;
    throw BRAINSCutExceptionStringHandler( errorMsg );
    }

  openCVANN->load( ANNModelFilename.c_str() );
}

inline scalarType *
BRAINSCutApplyModel
::GetArrayFromVector( scalarType array[], InputVectorType& vector, unsigned int inputVectorSize)
{
  for( unsigned int ai = 0; ai < inputVectorSize; ai++ )
    {
    array[ai] = vector[ai];
    }
  return array;
}

inline void
BRAINSCutApplyModel
::WritePredictROIProbabilityBasedOnReferenceImage( const PredictValueMapType& predictedOutput,
                                                   const WorkingImagePointer& referenceImage,
                                                   const WorkingImagePointer& roi,
                                                   const std::string imageFilename)
{
  WorkingImagePointer ANNContinuousOutputImage = WorkingImageType::New();

  ANNContinuousOutputImage->CopyInformation( referenceImage );
  ANNContinuousOutputImage->SetRegions( referenceImage->GetLargestPossibleRegion() );
  ANNContinuousOutputImage->Allocate();
  ANNContinuousOutputImage->FillBuffer( 0.0F );
  for( PredictValueMapType::const_iterator it = predictedOutput.begin();
       it != predictedOutput.end();
       ++it )
    {
    WorkingImageType::IndexType indexFromKey = FeatureInputVector::HashIndexFromKey( it->first );

    ANNContinuousOutputImage->SetPixel( indexFromKey, it->second );
    }

  itk::ImageRegionIterator<WorkingImageType> imgIt( roi, roi->GetLargestPossibleRegion() );
  imgIt.GoToBegin();
  while( !imgIt.IsAtEnd() )
    {
    if( imgIt.Value() > (HundreadPercentValue - FLOAT_TOLERANCE) )
      {
      ANNContinuousOutputImage->SetPixel( imgIt.GetIndex(), HundreadPercentValue );
      }
    ++imgIt;
    }

  itkUtil::WriteImage<WorkingImageType>( ANNContinuousOutputImage, imageFilename);
}

/* get output file dir */
inline std::string
BRAINSCutApplyModel
::GetSubjectOutputDirectory( DataSet& subject)
{
  std::string outputDir = subject.GetAttribute<StringValue>("OutputDir");

  if( !itksys::SystemTools::FileExists( outputDir.c_str(), false ) )
    {
    std::cout << " Subject output directory does not exist. Create as following."
              << outputDir.c_str()
              << std::endl;
    itksys::SystemTools::MakeDirectory( outputDir.c_str() );
    }
  return outputDir;
}

/* get continuous file name */
inline std::string
BRAINSCutApplyModel
::GetContinuousPredictionFilename( DataSet& subject, std::string currentROIName)
{
  const std::string subjectID(subject.GetAttribute<StringValue>("Name") );

  std::string outputDir = GetSubjectOutputDirectory( subject );

  std::string outputVolumeFilename = outputDir + "/ANNContinuousPrediction"
    + currentROIName
    + subjectID
    + ".nii.gz";

  return outputVolumeFilename;
}

/* get output mask file name of subject */
inline std::string
BRAINSCutApplyModel
::GetOutputROIFilename( DataSet& subject, std::string currentROIName)
{
  std::string givenROIName = subject.GetMaskFilenameByType( currentROIName );

  if( givenROIName == "" or givenROIName == "na" )
    {
    std::string outputDir =  GetSubjectOutputDirectory( subject );
    givenROIName = outputDir + "/ANNLabel_" + currentROIName + ".nii.gz";
    }
  return givenROIName;
}
