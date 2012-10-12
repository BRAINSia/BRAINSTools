#include "BRAINSCutApplyModel.h"
#include "FeatureInputVector.h"
#include "TrainingPrameters.h"
#include "ApplyModel.h"

#include <itkConnectedComponentImageFilter.h>
#include <itkRelabelComponentImageFilter.h>
#include <itkBinaryBallStructuringElement.h>
#include <itkBinaryMorphologicalClosingImageFilter.h>
#include <itkSigmoidImageFilter.h>

// TODO: consider using itk::LabelMap Hole filling process in ITK4

BRAINSCutApplyModel
::BRAINSCutApplyModel( BRAINSCutDataHandler& dataHandler )
  : numberOfTrees(-1),
  depthOfTree(-1)
{
  this->m_myDataHandler = &dataHandler;
  // TODO Take this apart to generate registration one by one!
  this->m_myDataHandler->SetRegionsOfInterest();
  this->m_myDataHandler->SetRegistrationParameters();

  // this->m_myDataHandler->SetAtlasDataSet();
  this->m_myDataHandler->SetRhoPhiTheta();

  this->m_myDataHandler->SetTrainingVectorConfiguration(); // this has to be before gradient size
  this->m_myDataHandler->SetGradientSize();

  this->m_myDataHandler->SetNormalization();

  gaussianSmoothingSigma = this->m_myDataHandler->GetGaussianSmoothingSigma();
  trainIteration   = this->m_myDataHandler->GetTrainIteration();
  applyDataSetList = this->m_myDataHandler->GetApplyDataSet();

  /** set default: ANN **/
  this->m_openCVANN = new OpenCVMLPType();

  SetMethod( "ANN" );
}

BRAINSCutApplyModel
::~BRAINSCutApplyModel()
{
  this->m_openCVANN->clear();
  delete this->m_openCVANN;

  if(  method ==  "RandomForest" )
    {
    this->openCVRandomForest->clear();
    delete this->openCVRandomForest;
    }
}

void
BRAINSCutApplyModel
::SetMethod( std::string inputMethod)
{
  method = inputMethod;
  if( method ==  "ANN" )
    {
    this->m_myDataHandler->SetTrainConfiguration( "ANNParameters" );
    }
  else
    {
    openCVRandomForest = new CvRTrees;
    this->m_myDataHandler->SetTrainConfiguration( "RandomForestParameters");
    }
}

/* iterate through subject */
void
BRAINSCutApplyModel
::Apply()
{
  if( method ==  "ANN" )
    {
    this->m_myDataHandler->SetANNTestingSSEFilename();
    annOutputThreshold = this->m_myDataHandler->GetANNOutputThreshold();
    this->m_myDataHandler->SetANNModelFilenameAtIteration( trainIteration);
    ReadANNModelFile();
    }
  else if( method == "RandomForest" )
    {
    if( this->m_myDataHandler->GetRandomForestModelFilename() == "" )
      {
      this->m_myDataHandler->SetRandomForestModelFilename( depthOfTree, numberOfTrees);
      }
    ReadRandomForestModelFile();
    }

  typedef BRAINSCutConfiguration::ApplyDataSetListType::iterator ApplySubjectIteratorType;

  normalization = this->m_myDataHandler->GetNormalization();
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

  this->m_myDataHandler->GetDeformedSpatialLocationImages( deformedSpatialLocationImageList, subject );

  WorkingImageVectorType imagesOfInterest;
  this->m_myDataHandler->ReadImagesOfSubjectInOrder(imagesOfInterest, subject);

  /** Warp probability map(ROI) onto the subject*/
  typedef std::map<std::string, WorkingImagePointer> DeformedROIMapType;
  DeformedROIMapType deformedROIs;

  this->m_myDataHandler->GetDeformedROIs(deformedROIs, subject);

  /** Gaussian Smoothing if requested to cover broader area */
  if( gaussianSmoothingSigma > 0.0F )
    {
    for( DeformedROIMapType::iterator it = deformedROIs.begin();
         it != deformedROIs.end();
         ++it )
      {
      deformedROIs[it->first] = SmoothImage( it->second, gaussianSmoothingSigma);
      }
    }

  /** Get input feature vectors based on subject images and deformed ROIs */
  FeatureInputVector inputVectorGenerator;

  inputVectorGenerator.SetGradientSize( this->m_myDataHandler->GetGradientSize() );
  inputVectorGenerator.SetImagesOfInterestInOrder( imagesOfInterest );
  inputVectorGenerator.SetImagesOfSpatialLocation( deformedSpatialLocationImageList );
  inputVectorGenerator.SetCandidateROIs( deformedROIs);
  inputVectorGenerator.SetROIInOrder( this->m_myDataHandler->GetROIIDsInOrder() );
  inputVectorGenerator.SetInputVectorSize();
  inputVectorGenerator.SetNormalization( normalization );

  /* now iterate through the roi */

  unsigned int roiIDsOrderNumber = 0;
  // for( DataSet::StringVectorType::iterator roiTyIt = this->m_myDataHandler->GetROIIDsInOrder().begin();
  //     roiTyIt != this->m_myDataHandler->GetROIIDsInOrder().end();
  //     ++roiTyIt ) // roiTyIt = Region of Interest Type Iterator
  while( roiIDsOrderNumber < this->m_myDataHandler->GetROIIDsInOrder().size() )
    {
    const std::string currentROIName = std::string( this->m_myDataHandler->GetROIIDsInOrder()[roiIDsOrderNumber] );

    ProbabilityMapParser* roiDataSet =
      this->m_myDataHandler->GetROIDataList()->GetMatching<ProbabilityMapParser>( "StructureID",
                                                                                  currentROIName.c_str() );
    if( roiDataSet->GetAttribute<StringValue>("GenerateVector") == "true" )
      {
      InputVectorMapType  roiInputVector = inputVectorGenerator.ComputeAndGetFeatureInputOfROI( currentROIName );
      PredictValueMapType predictedOutputVector;

      if( !computeSSE )
        {
        PredictROI( roiInputVector, predictedOutputVector,
                    roiIDsOrderNumber, inputVectorGenerator.GetInputVectorSize() );
        roiInputVector.clear();
        const std::string ANNContinuousOutputFilename = GetContinuousPredictionFilename( subject, currentROIName );

        /* post processing
         * may include hole-filling(closing), thresholding, and more adjustment
         */
        BinaryImagePointer mask;
        if( method == "ANN" )
          {
          WritePredictROIProbabilityBasedOnReferenceImage( predictedOutputVector,
                                                           imagesOfInterest.front(),
                                                           deformedROIs.find( currentROIName )->second,
                                                           ANNContinuousOutputFilename,
                                                           1.0F );
          mask = PostProcessingANN( ANNContinuousOutputFilename,
                                    annOutputThreshold);
          }
        else if( method == "RandomForest" )
          {
          WritePredictROIProbabilityBasedOnReferenceImage( predictedOutputVector,
                                                           imagesOfInterest.front(),
                                                           deformedROIs.find( currentROIName )->second,
                                                           ANNContinuousOutputFilename,
                                                           roiIDsOrderNumber + 1 );
          mask = PostProcessingRF( ANNContinuousOutputFilename );
          }

        std::string roiOutputFilename = GetROIVolumeName( subject, currentROIName );
        itkUtil::WriteImage<BinaryImageType>( mask, roiOutputFilename );
        itkUtil::WriteImage<WorkingImageType>( deformedROIs.find( currentROIName )->second,
                                               roiOutputFilename + "def.nii.gz");
        }
      else /* testing phase */
        {
        for( int currentIteration = 1; currentIteration <= trainIteration; currentIteration++ )
          {
          this->m_myDataHandler->SetANNModelFilenameAtIteration( currentIteration );
          PredictROI( roiInputVector, predictedOutputVector,
                      roiIDsOrderNumber, inputVectorGenerator.GetInputVectorSize() );
          roiInputVector.clear();
          const std::string roiReferenceFilename = GetROIVolumeName( subject, currentROIName );
          const float       SSE = ComputeSSE( predictedOutputVector, roiReferenceFilename );

          ANNTestingSSEFileStream << currentROIName
                                  << ", subjectID, " << subjectID
                                  << ", Iteration, " << currentIteration
                                  << ", SSE, " << SSE
                                  << std::endl;
          }
        }

      predictedOutputVector.clear();
      }
    roiIDsOrderNumber++;
    }

  deformedSpatialLocationImageList.clear();
  imagesOfInterest.clear();
  deformedROIs.clear();
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

BinaryImagePointer
BRAINSCutApplyModel
::PostProcessingANN( std::string continuousFilename,
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
::PostProcessingRF( std::string labelImageFilename )
{
  WorkingImagePointer labelImage = ReadImageByFilename( labelImageFilename );

  BinaryImagePointer maskVolume = itkUtil::ScaleAndCast<WorkingImageType,
                                                        BinaryImageType>( labelImage,
                                                                          0,
                                                                          255);

  /* Get One label */
  BinaryImagePointer connected_maskVolume = GetOneConnectedRegion( maskVolume );

  /* opening and closing to get rid of island and holes */
  BinaryImagePointer filled_maskVolume = FillHole( connected_maskVolume );

  return filled_maskVolume;
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

  BinaryImagePointer mask = itkUtil::ScaleAndCast<WorkingImageType,
                                                  BinaryImageType>(thresholder->GetOutput(),
                                                                   0,
                                                                   255);
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
  thresholder->SetUpperThreshold( 255 );
  thresholder->Update();

  BinaryImagePointer mask = itkUtil::ScaleAndCast<WorkingImageType,
                                                  BinaryImageType>(thresholder->GetOutput(),
                                                                   0,
                                                                   255);
  return mask;
}

BinaryImagePointer
BRAINSCutApplyModel
::ExtractLabel( BinaryImagePointer& image, unsigned char thresholdValue  )
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
::GetOneConnectedRegion( BinaryImagePointer& image )
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

BinaryImagePointer
BRAINSCutApplyModel
::FillHole( BinaryImagePointer& mask)
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
              const unsigned int         roiNumber,
              const unsigned int         inputVectorSize) const
{
  const unsigned int NumberOfROIS = this->m_myDataHandler->GetROIIDsInOrder().size();

  // The openCVOutput is unnecessary of RandomForest, but it is really cheap to do once.
  // even if it is wasted effort.
  // CvMat * openCVOutput = cvCreateMat( 1, NumberOfROIS, CV_32FC1);

  /* initialize container of output vector*/
  resultOutputVector.clear();
  for( InputVectorMapType::iterator it = roiInputFeatureVector.begin();
       it != roiInputFeatureVector.end();
       ++it )
    {
    /* get open cv type matrix from array for prediction */

    CvMat openCVInputFeature;
    cvInitMatHeader( &openCVInputFeature, 1, inputVectorSize, CV_32FC1,
                     &(it->second[0]) ); //
                                         // http://stackoverflow.com/questions/6701816/how-to-use-a-stdvector-in-a-c-function
    // std::cout<<FeatureInputVector::HashIndexFromKey( it->first );

    /* predict */
    if( method == "ANN" )
      {
      CvMat * openCVOutput = cvCreateMat( 1, NumberOfROIS, CV_32FC1 );
      this->m_openCVANN->predict( &openCVInputFeature, openCVOutput );

      /* insert result to the result output vector */
      resultOutputVector.insert( std::pair<hashKeyType, scalarType>(
                                   ( it->first ),
                                   CV_MAT_ELEM( *openCVOutput, scalarType, 0, roiNumber) ) );
      }
    else if( method == "RandomForest" )
      {
      const scalarType response = openCVRandomForest->predict( &openCVInputFeature );
      // make binary input
      /*if( response > 0.5F )
      {
        response = 1.0F;
      }*/
      resultOutputVector.insert( std::pair<hashKeyType, scalarType>(  ( it->first ),
                                                                      response ) );
      // std::cout<<"--> "<<response<<std::endl;
      }
    }
}

void
BRAINSCutApplyModel
::SetComputeSSE( const bool sse)
{
  computeSSE = sse;
  if( computeSSE )
    {
    ANNTestingSSEFileStream.open( this->m_myDataHandler->GetANNTestingSSEFilename().c_str(),
                                  std::fstream::out );
    if( !ANNTestingSSEFileStream.good() )
      {
      std::string errorMsg = " Cannot open the file! :";
      errorMsg += this->m_myDataHandler->GetANNTestingSSEFilename();
      std::cout << errorMsg << std::endl;
      throw BRAINSCutExceptionStringHandler( errorMsg );
      exit( EXIT_FAILURE );
      }
    }
}

/** read model files */
void
BRAINSCutApplyModel
::ReadANNModelFile()
{
  std::string ANNModelFilename = this->m_myDataHandler->GetANNModelFilename();

  if( !itksys::SystemTools::FileExists( ANNModelFilename.c_str() ) )
    {
    std::string errorMsg = " File does not exist! :";
    errorMsg += ANNModelFilename;
    throw BRAINSCutExceptionStringHandler( errorMsg );
    }

  std::cout << "Filename:: " << ANNModelFilename << std::endl;
  this->m_openCVANN->load( ANNModelFilename.c_str() );
}

void
BRAINSCutApplyModel
::ReadRandomForestModelFile()
{
  const std::string randomForestFilename = this->m_myDataHandler->GetRandomForestModelFilename();

  if( !itksys::SystemTools::FileExists( randomForestFilename.c_str() ) )
    {
    std::string errorMsg = " File does not exist! :";
    errorMsg += randomForestFilename;
    throw BRAINSCutExceptionStringHandler( errorMsg );
    }
  std::cout << "******* LOAD random forest file ********" << std::endl;
  openCVRandomForest->load( randomForestFilename.c_str() );
}

inline void
BRAINSCutApplyModel
::WritePredictROIProbabilityBasedOnReferenceImage( const PredictValueMapType& predictedOutput,
                                                   const WorkingImagePointer& referenceImage,
                                                   const WorkingImagePointer& roi,
                                                   const std::string imageFilename,
                                                   const WorkingPixelType labelValue )
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
    if( imgIt.Value() >= (HundredPercentValue - FLOAT_TOLERANCE) )
      {
      ANNContinuousOutputImage->SetPixel( imgIt.GetIndex(), labelValue );
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

void
BRAINSCutApplyModel
::SetNumberOfTrees( const int trees)
{
  if( trees < 0 )
    {
    numberOfTrees = this->m_myDataHandler->GetMaxTreeCount();
    }
  else
    {
    numberOfTrees = trees;
    }
}

void
BRAINSCutApplyModel
::SetDepthOfTree( const int depth )
{
  if( depth < 0 )
    {
    depthOfTree = this->m_myDataHandler->GetMaxDepth();
    }
  else
    {
    depthOfTree = depth;
    }
}
