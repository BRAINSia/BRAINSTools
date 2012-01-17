#include "BRAINSCutApplyModel.h"
#include "FeatureInputVector.h"
#include "ANNParams.h"
#include "ApplyModel.h"

#include <itkConnectedComponentImageFilter.h>
#include <itkRelabelComponentImageFilter.h>
#include <itkBinaryBallStructuringElement.h>
#include <itkBinaryMorphologicalClosingImageFilter.h>
#include <itkSigmoidImageFilter.h>

// TODO: consider using itk::LabelMap Hole filling process in ITK4

BRAINSCutApplyModel
::BRAINSCutApplyModel( std::string netConfigurationFilename)
  :
  BRAINSCutPrimary( netConfigurationFilename )
  // ANNModelFilename(NULL)
{
  // TODO Take this apart to generate registration one by one!
  SetRegionsOfInterestFromNetConfiguration();
  SetRegistrationParametersFromNetConfiguration();

  SetAtlasDataSet();
  SetRhoPhiThetaFromNetConfiguration();

  SetANNModelConfiguration();
  SetANNTestingSSEFilename();
  SetGradientSizeFromNetConfiguration();
  SetANNOutputThresholdFromNetConfiguration();
  SetTrainIterationFromNetConfiguration();

  SetApplyDataSetFromNetConfiguration();
  SetANNModelFilenameFromNetConfiguration();
  SetGaussianSmoothingSigmaFromNetConfiguration();

  openCVANN = new OpenCVMLPType();
}

/* iterate through subject */
void
BRAINSCutApplyModel
::Apply()
{
  typedef NetConfiguration::ApplyDataSetListType::iterator ApplySubjectIteratorType;

  normalization = GetNormalizationFromNetConfiguration();
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
  const std::string subjectID(subject.GetAttribute<StringValue>("Name") );

  std::map<std::string, WorkingImagePointer> deformedSpatialLocationImageList;

  GetDeformedSpatialLocationImages( deformedSpatialLocationImageList, subject );

  WorkingImageVectorType imagesOfInterest;
  GetImagesOfSubjectInOrder(imagesOfInterest, subject);

  /** Warp probability map(ROI) onto the subject*/
  typedef std::map<std::string, WorkingImagePointer> DeformedROIMapType;
  DeformedROIMapType deformedROIs;

  GetDeformedROIs(deformedROIs, subject);

  /** Gaussian Smoothing if requested to cover broader area */
  if( gaussianSmoothingSigma > 0.0F )
    {
    for( DeformedROIMapType::iterator it = deformedROIs.begin();
         it != deformedROIs.end();
         ++it )
      {
      deformedROIs[it->first] = this->SmoothImage( it->second, gaussianSmoothingSigma);
      }
    }

  /** Get input feature vectors based on subject images and deformed ROIs */
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

      if( !computeSSE )
        {
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
        BinaryImagePointer mask = PostProcessingOfANNContinuousImage( ANNContinuousOutputFilename,
                                                                      annOutputThreshold);

        std::string roiOutputFilename = GetROIVolumeName( subject, *roiTyIt );
        itkUtil::WriteImage<BinaryImageType>( mask, roiOutputFilename );
        itkUtil::WriteImage<WorkingImageType>( deformedROIs.find( *roiTyIt )->second,
                                               roiOutputFilename + "def.nii.gz");
        }
      else /* testing phase */
        {
        for( int currentIteration = 1; currentIteration <= trainIteration; currentIteration++ )
          {
          SetANNModelFilenameAtIteration( currentIteration );
          PredictROI( roiInputVector, predictedOutputVector,
                      roiIDsOrderNumber, inputVectorGenerator.GetInputVectorSize() );
          std::string roiReferenceFilename = GetROIVolumeName( subject, *roiTyIt );
          float       SSE = ComputeSSE( predictedOutputVector, roiReferenceFilename );

          ANNTestingSSEFileStream << *roiTyIt
                                  << ", subjectID, " << subjectID
                                  << ", Iteration, " << currentIteration
                                  << ", SSE, " << SSE
                                  << std::endl;
          }
        }
      }

    roiIDsOrderNumber++;
    }
}

float
BRAINSCutApplyModel
::ComputeSSE( const PredictValueMapType& predictedOutputVector,
              const std::string roiReferenceFilename )
{
  WorkingImagePointer ReferenceVolume = ReadImageByFilename( roiReferenceFilename );

  WorkingImageType::PixelType referenceValue = 0.0F;
  double                      SSE = 0.0F;

  for( PredictValueMapType::const_iterator it = predictedOutputVector.begin();
       it != predictedOutputVector.end();
       ++it )
    {
    WorkingImageType::IndexType indexFromKey = FeatureInputVector::HashIndexFromKey( it->first );
    referenceValue = ReferenceVolume->GetPixel( indexFromKey );
    SSE += (referenceValue - it->second) * (referenceValue - it->second);
    }
  double totalSize = predictedOutputVector.size();
  SSE = SSE / totalSize;
  return SSE;
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

BinaryImagePointer
BRAINSCutApplyModel
::PostProcessingOfANNContinuousImage( std::string continuousFilename,
                                      scalarType threshold )
{
  WorkingImagePointer continuousImage = ReadImageByFilename( continuousFilename );

  BinaryImagePointer maskVolume;

  /* threshold */
  maskVolume = ThresholdImageAtLower( continuousImage, threshold);

  /* Get One label */
  maskVolume = GetOneConnectedRegion( maskVolume );

  /* opening and closing to get rid of island and holes */
  maskVolume = FillHole( maskVolume );
  return maskVolume;
}

BinaryImagePointer
BRAINSCutApplyModel
::ThresholdImageAtUpper( WorkingImagePointer& image, scalarType thresholdValue  )
{
  typedef itk::BinaryThresholdImageFilter<WorkingImageType, WorkingImageType> ThresholdFilterType;
  ThresholdFilterType::Pointer thresholder = ThresholdFilterType::New();

  if( thresholdValue < 0.0F )
    {
    std::string msg = " ANNOutput Threshold cannot be less than zero. \n";
    throw BRAINSCutExceptionStringHandler( msg );
    }
  thresholder->SetInput( image );
  thresholder->SetOutsideValue( 0 );
  thresholder->SetInsideValue( 1 );
  thresholder->SetUpperThreshold( thresholdValue );
  thresholder->SetLowerThreshold( -1e+10F);
  thresholder->Update();

  BinaryImagePointer mask = itkUtil::TypeCast<WorkingImageType, BinaryImageType>( thresholder->GetOutput() );
  return mask;
}

BinaryImagePointer
BRAINSCutApplyModel
::ThresholdImageAtLower( WorkingImagePointer& image, scalarType thresholdValue  )
{
  typedef itk::BinaryThresholdImageFilter<WorkingImageType, WorkingImageType> ThresholdFilterType;
  ThresholdFilterType::Pointer thresholder = ThresholdFilterType::New();

  if( thresholdValue < 0.0F )
    {
    std::string msg = " ANNOutput Threshold cannot be less than zero. \n";
    throw BRAINSCutExceptionStringHandler( msg );
    }
  thresholder->SetInput( image );
  thresholder->SetOutsideValue( 0 );
  thresholder->SetInsideValue( 1 );
  thresholder->SetLowerThreshold( thresholdValue );
  thresholder->Update();

  BinaryImagePointer mask = itkUtil::TypeCast<WorkingImageType, BinaryImageType>( thresholder->GetOutput() );
  return mask;
}

BinaryImagePointer
BRAINSCutApplyModel
::ExtractLabel( BinaryImagePointer image, unsigned char thresholdValue  )
{
  typedef itk::BinaryThresholdImageFilter<BinaryImageType, BinaryImageType> ThresholdFilterType;
  ThresholdFilterType::Pointer thresholder = ThresholdFilterType::New();

  thresholder->SetInput( image );
  thresholder->SetInsideValue( 1 );
  thresholder->SetOutsideValue( 0 );
  thresholder->SetUpperThreshold( thresholdValue );
  thresholder->SetLowerThreshold( thresholdValue );
  thresholder->Update();

  BinaryImagePointer mask = itkUtil::TypeCast<BinaryImageType, BinaryImageType>( thresholder->GetOutput() );
  return mask;
}

BinaryImagePointer
BRAINSCutApplyModel
::GetOneConnectedRegion( BinaryImagePointer image )
{
  /* relabel images if they are disconnected */
  typedef itk::ConnectedComponentImageFilter<BinaryImageType, BinaryImageType>
    ConnectedBinaryImageFilterType;
  ConnectedBinaryImageFilterType::Pointer relabler = ConnectedBinaryImageFilterType::New();

  relabler->SetInput( image );

  /* relable images from the largest to smallset one */
  typedef itk::RelabelComponentImageFilter<BinaryImageType, BinaryImageType> RelabelInOrderFilterType;
  RelabelInOrderFilterType::Pointer relabelInOrder = RelabelInOrderFilterType::New();

  relabelInOrder->SetInput( relabler->GetOutput() );
  relabelInOrder->Update();

  BinaryImagePointer multipleLabelVolume = relabelInOrder->GetOutput();
  /* get the label one */
  BinaryImagePointer resultMask = ExtractLabel( multipleLabelVolume, 1 );

  return resultMask;
}

void
BRAINSCutApplyModel
::SetGaussianSmoothingSigmaFromNetConfiguration()
{
  gaussianSmoothingSigma =
    BRAINSCutNetConfiguration.Get<ApplyModelType>("ApplyModel")->GetAttribute<FloatValue>("GaussianSmoothingSigma");
}

void
BRAINSCutApplyModel
::SetANNOutputThresholdFromNetConfiguration()
{
  annOutputThreshold =
    BRAINSCutNetConfiguration.Get<ApplyModelType>("ApplyModel")->GetAttribute<FloatValue>("MaskThresh");
  if( annOutputThreshold < 0.0F )
    {
    std::string msg = " ANNOutput Threshold cannot be less than zero. \n";
    throw BRAINSCutExceptionStringHandler( msg );
    }
}

BinaryImagePointer
BRAINSCutApplyModel
::FillHole( BinaryImagePointer mask)
{
  typedef itk::BinaryBallStructuringElement<BinaryImageType::PixelType, DIMENSION> KernelType;
  KernelType           ball;
  KernelType::SizeType ballSize;

  /* Create the structuring element- a disk of radius 2 */
  ballSize.Fill(2);
  ball.SetRadius( ballSize );
  ball.CreateStructuringElement();

  /* Closing */
  typedef itk::BinaryMorphologicalClosingImageFilter<BinaryImageType,
                                                     BinaryImageType,
                                                     KernelType> ClosingFilterType;
  ClosingFilterType::Pointer closingFilter = ClosingFilterType::New();
  closingFilter->SetInput( mask );
  closingFilter->SetKernel( ball );
  closingFilter->SetForegroundValue(1);
  closingFilter->SetSafeBorder( true );
  closingFilter->Update();

  return closingFilter->GetOutput();
}

inline void
BRAINSCutApplyModel
::PredictROI( InputVectorMapType&  roiInputFeatureVector,
              PredictValueMapType& resultOutputVector,
              unsigned int         roiNumber,
              unsigned int         inputVectorSize)
{
  ReadANNModelFile();
  /* initialize container of output vector*/
  resultOutputVector.clear();
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
    resultOutputVector.insert( std::pair<int, scalarType>(  ( it->first ),  CV_MAT_ELEM( *openCVOutput,
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

std::string
BRAINSCutApplyModel
::GetANNModelBaseName()
{
  std::string basename;

  try
    {
    basename = annModelConfiguration->GetAttribute<StringValue>("TrainingModelFilename");
    }
  catch( ... )
    {
    throw BRAINSCutExceptionStringHandler("Fail to get the ann model file name");
    }
  return basename;
}

void
BRAINSCutApplyModel
::SetANNModelFilenameAtIteration( const int iteration)
{
  ANNModelFilename = GetANNModelBaseName();

  char temp[10];
  sprintf( temp, "%09d", iteration );
  ANNModelFilename += temp;
}

void
BRAINSCutApplyModel
::SetTrainIterationFromNetConfiguration()
{
  trainIteration = BRAINSCutNetConfiguration.Get<ANNParams>("ANNParams")->GetAttribute<IntValue>("Iterations");
}

void
BRAINSCutApplyModel
::SetANNTestingSSEFilename()
{
  ANNTestingSSEFilename = GetANNModelBaseName();
  ANNTestingSSEFilename += "ValidationSetSSE.txt";
}

void
BRAINSCutApplyModel
::SetComputeSSE( const bool sse)
{
  computeSSE = sse;
  if( computeSSE )
    {
    ANNTestingSSEFileStream.open( ANNTestingSSEFilename.c_str(),
                                  std::fstream::out );
    if( !ANNTestingSSEFileStream.good() )
      {
      std::string errorMsg = " Cannot open the file! :";
      errorMsg += ANNTestingSSEFilename;
      std::cout << errorMsg << std::endl;
      throw BRAINSCutExceptionStringHandler( errorMsg );
      exit( EXIT_FAILURE );
      }
    }
}

void
BRAINSCutApplyModel
::SetANNModelFilenameFromNetConfiguration()
{
  SetANNModelFilenameAtIteration( trainIteration );
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
::GetROIVolumeName( DataSet& subject, std::string currentROIName)
{
  std::string       givenROIName = subject.GetMaskFilenameByType( currentROIName );
  const std::string subjectID(subject.GetAttribute<StringValue>("Name") );

  if( givenROIName == "" or givenROIName == "na" )
    {
    std::string outputDir =  GetSubjectOutputDirectory( subject );
    givenROIName = outputDir + "/" + subjectID + "ANNLabel_" + currentROIName + ".nii.gz";
    }
  return givenROIName;
}
