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
#include "BRAINSCutApplyModel.h"
#include "BRAINSCutUtilities.h"
#include "FeatureInputVector.h"
#include "TrainingPrameters.h"
#include "ApplyModel.h"

#include <itkConnectedComponentImageFilter.h>
#include <itkRelabelComponentImageFilter.h>
#include <itkBinaryBallStructuringElement.h>
#include <itkBinaryMorphologicalClosingImageFilter.h>
#include <itkBinaryMorphologicalOpeningImageFilter.h>
#include <itkSigmoidImageFilter.h>
#include <itkImageDuplicator.h>
#include <itkLabelStatisticsImageFilter.h>
#include <itkBinaryFillholeImageFilter.h>
#include <itkMultiplyImageFilter.h>
#include "itkNumberToString.h"

// TODO: consider using itk::LabelMap Hole filling process in ITK4
BRAINSCutApplyModel
::BRAINSCutApplyModel() :
  m_myDataHandler(nullptr),
  m_applyDataSetList(),
  m_method(""),
  m_normalization(""),
  m_computeSSE(true),
  m_trainIteration(0),
  m_numberOfTrees(0),
  m_depthOfTree(0),
  m_ANNTestingSSEFileStream(),
  m_annOutputThreshold(0),
  m_gaussianSmoothingSigma(0)
{
}

#if 0
BRAINSCutApplyModel
::BRAINSCutApplyModel( BRAINSCutApplyModel& applyModel ) :
  m_myDataHandler(NULL),
  m_applyDataSetList(),
  m_method(""),
  m_normalization(""),
  m_computeSSE(true),
  m_trainIteration(0),
  m_numberOfTrees(0),
  m_depthOfTree(0),
  m_ANNTestingSSEFileStream(),
  m_annOutputThreshold(0),
  m_gaussianSmoothingSigma(0),
  m_openCVANN(NULL),
  m_openCVRandomForest(NULL)
{
}
#endif

BRAINSCutApplyModel
::BRAINSCutApplyModel( BRAINSCutDataHandler& dataHandler ) :
  m_myDataHandler(nullptr),
  m_applyDataSetList(),
  m_method(""),
  m_normalization(""),
  m_computeSSE(true),
  m_trainIteration(0),
  m_numberOfTrees(0),
  m_depthOfTree(0),
  m_ANNTestingSSEFileStream(),
  m_annOutputThreshold(0),
  m_gaussianSmoothingSigma(0)
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

  m_gaussianSmoothingSigma = this->m_myDataHandler->GetGaussianSmoothingSigma();
  m_trainIteration   = this->m_myDataHandler->GetTrainIteration();
  m_applyDataSetList = this->m_myDataHandler->GetApplyDataSet();

  /** set default: ANN **/
  this->m_openCVANN = OpenCVMLPType::create();
  //new OpenCVMLPType();

  SetMethod( "ANN" );
}

BRAINSCutApplyModel
::~BRAINSCutApplyModel()
{
}

void
BRAINSCutApplyModel
::SetMethod( std::string inputMethod)
{
  m_method = inputMethod;
  if( m_method ==  "ANN" )
    {
    this->m_myDataHandler->SetTrainConfiguration( "ANNParameters" );
    }
  else
    {
    //m_openCVRandomForest = CvRTrees;
    m_openCVRandomForest = cv::ml::RTrees::create();
    this->m_myDataHandler->SetTrainConfiguration( "RandomForestParameters");
    }
}

/* iterate through subject */
void
BRAINSCutApplyModel
::Apply()
{
  if( m_method ==  "ANN" )
    {
    this->m_myDataHandler->SetANNTestingSSEFilename();
    m_annOutputThreshold = this->m_myDataHandler->GetANNOutputThreshold();
    this->m_myDataHandler->SetANNModelFilenameAtIteration( m_trainIteration);
    ReadANNModelFile();
    }
  else if( m_method == "RandomForest" )
    {
    if( this->m_myDataHandler->GetRandomForestModelFilename() == "" )
      {
      this->m_myDataHandler->SetRandomForestModelFilename( m_depthOfTree, m_numberOfTrees);
      }
    ReadRandomForestModelFile();
    }

  typedef BRAINSCutConfiguration::ApplyDataSetListType::iterator ApplySubjectIteratorType;

  m_normalization = this->m_myDataHandler->GetNormalizationMethod();
  for( ApplySubjectIteratorType subjectIt = m_applyDataSetList.begin();
       subjectIt != m_applyDataSetList.end();
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
  itk::NumberToString<double>    doubleToString;
  const std::string subjectID(subject.GetAttribute<StringValue>("Name") );

  std::map<std::string, WorkingImagePointer> deformedSpatialLocationImageList;

  this->m_myDataHandler->GetDeformedSpatialLocationImages( deformedSpatialLocationImageList, subject );

  WorkingImageVectorType imagesOfInterest;
  this->m_myDataHandler->ReadImagesOfSubjectInOrder(imagesOfInterest, subject);

  /** Warp probability map(ROI) onto the subject*/
  typedef std::map<std::string, WorkingImagePointer> DeformedROIMapType;
  DeformedROIMapType deformedROIs;

  this->m_myDataHandler->GetDeformedROIs(deformedROIs, subject);

  WorkingImagePointer candidateRegion =  m_myDataHandler->GetCandidateRegion( subject );

  if( candidateRegion.IsNotNull() )
    {
    // TODO Clip deformed ROI here
    for( DeformedROIMapType::iterator it = deformedROIs.begin();
         it != deformedROIs.end();
         ++it )
      {
      it->second = ClipImageWithBinaryMask( it->second, candidateRegion );
      }
    }

  /** Gaussian Smoothing if requested to cover broader area */
  if( m_gaussianSmoothingSigma > 0.0F )
    {
    for( DeformedROIMapType::iterator it = deformedROIs.begin();
         it != deformedROIs.end();
         ++it )
      {
      deformedROIs[it->first] = SmoothImage( it->second, m_gaussianSmoothingSigma);
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
  inputVectorGenerator.SetNormalizationMethod( m_normalization );

  /* now iterate through the roi */

  unsigned int          roiIDsOrderNumber = 0;
  LabelImagePointerType resultLabelFromRF;
  LabelImagePointerType ambiguousLabelFromRF;

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

      if( !m_computeSSE )
        {
        PredictROI( roiInputVector, predictedOutputVector,
                    roiIDsOrderNumber, inputVectorGenerator.GetInputVectorSize() );
        roiInputVector.clear();
        const std::string & ANNContinuousOutputFilename = GetContinuousPredictionFilename( subject, currentROIName );

        /* post processing
         * may include hole-filling(closing), thresholding, and more adjustment
         */

        LabelImagePointerType mask;
        std::string           roiOutputFilename = GetROIVolumeName( subject, currentROIName );
        if( m_method == "ANN" )
          {
          WritePredictROIProbabilityBasedOnReferenceImage( predictedOutputVector,
                                                           imagesOfInterest.front(),
                                                           deformedROIs.find( currentROIName )->second,
                                                           ANNContinuousOutputFilename,
                                                           1.0F );
          mask = PostProcessingANN( ANNContinuousOutputFilename, this->m_annOutputThreshold);

          try
            {
            itkUtil::WriteImage<LabelImageType>( mask, roiOutputFilename );
            }
          catch( itk::ExceptionObject& ex )
            {
            std::cout << ex.what() << std::endl;
            }
          }
        else if( m_method == "RandomForest" )
          {
          WritePredictROIProbabilityBasedOnReferenceImage( predictedOutputVector,
                                                           imagesOfInterest.front(),
                                                           deformedROIs.find( currentROIName )->second,
                                                           ANNContinuousOutputFilename,
                                                           roiIDsOrderNumber + 1 );

          mask = PostProcessingRF( ANNContinuousOutputFilename, roiIDsOrderNumber + 1 );

          if( roiIDsOrderNumber == 0 )
            {
            typedef itk::ImageDuplicator<LabelImageType> DuplicatorType;
            DuplicatorType::Pointer labelDuplicator = DuplicatorType::New();
            labelDuplicator->SetInputImage( mask );
            labelDuplicator->Update();
            resultLabelFromRF = labelDuplicator->GetModifiableOutput();

            ambiguousLabelFromRF = LabelImageType::New();
            ambiguousLabelFromRF->CopyInformation( mask );
            ambiguousLabelFromRF->SetRegions( mask->GetLargestPossibleRegion() );
            ambiguousLabelFromRF->Allocate();
            ambiguousLabelFromRF->FillBuffer( 0 );
            }
          else
            {
            resultLabelFromRF = CombineLabel( resultLabelFromRF, mask, roiIDsOrderNumber + 1 );
            }
          ambiguousLabelFromRF = AmbiguousCountLabel( ambiguousLabelFromRF,
                                                      resultLabelFromRF,
                                                      mask);
          }

        itkUtil::WriteImage<WorkingImageType>( deformedROIs.find( currentROIName )->second,
                                               roiOutputFilename + "def.nii.gz");
        }
      else /* testing phase */
        {
        for( int currentIteration = 1; currentIteration <= m_trainIteration; ++currentIteration )
          {
          this->m_myDataHandler->SetANNModelFilenameAtIteration( currentIteration );
          PredictROI( roiInputVector, predictedOutputVector,
                      roiIDsOrderNumber, inputVectorGenerator.GetInputVectorSize() );
          roiInputVector.clear();
          const std::string roiReferenceFilename = GetROIVolumeName( subject, currentROIName );
          const float       SSE = ComputeSSE( predictedOutputVector, roiReferenceFilename );

          m_ANNTestingSSEFileStream << currentROIName
                                    << ", subjectID, " << subjectID
                                    << ", Iteration, " << currentIteration
                                    << ", SSE, " << doubleToString(SSE)
                                    << std::endl;
          }
        }
      predictedOutputVector.clear();
      }
    ++roiIDsOrderNumber;
    }

  if( m_method == "RandomForest" )
    {
    std::string labelFilename = GetLabelMapFilename( subject );
    WriteLabelMapToBinaryImages( subject, resultLabelFromRF );

    itkUtil::WriteImage<LabelImageType>( resultLabelFromRF, labelFilename );
    itkUtil::WriteImage<LabelImageType>( ambiguousLabelFromRF, labelFilename + "_AmbiguousMap.nii.gz");
    }

  deformedSpatialLocationImageList.clear();
  imagesOfInterest.clear();
  deformedROIs.clear();
}

float
BRAINSCutApplyModel
::ComputeSSE( const PredictValueMapType& predictedOutputVector,
              const std::string & roiReferenceFilename )
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
::WriteLabelMapToBinaryImages( const DataSet& subject,
                               const LabelImagePointerType& labelMapImage )
{
  try
    {
    for( unsigned char roiID = 0;
         roiID < this->m_myDataHandler->GetROIIDsInOrder().size();
         ++roiID )
      {
      LabelImagePointerType currentBinaryImage = ExtractLabel( labelMapImage,
                                                               roiID + 1  ); // label starts from 1
      std::string currentROIName = std::string( this->m_myDataHandler->GetROIIDsInOrder()[roiID] );
      std::string roiOutputFilename = GetROIVolumeName( subject, currentROIName );
      itkUtil::WriteImage<LabelImageType>( currentBinaryImage,
                                           roiOutputFilename );
      }
    }
  catch( std::exception& ex )
    {
    std::cout << ex.what() << std::endl;
    std::exit( EXIT_FAILURE );
    }
  return;
}

LabelImagePointerType
BRAINSCutApplyModel
::AmbiguousCountLabel( LabelImagePointerType& ambiguousMap,
                       LabelImagePointerType& combinedLabel,
                       LabelImagePointerType& currentLabel )
{
  try
    {
    typedef itk::ImageRegionIterator<LabelImageType> RegionIteratorType;

    RegionIteratorType itOut( ambiguousMap, ambiguousMap->GetLargestPossibleRegion() );
    RegionIteratorType itIn( currentLabel, currentLabel->GetLargestPossibleRegion() );
    RegionIteratorType itCompare( combinedLabel, combinedLabel->GetLargestPossibleRegion() );

    itIn.GoToBegin();
    itOut.GoToBegin();
    itCompare.GoToBegin();
    for( ; !itIn.IsAtEnd(); ++itCompare, ++itIn, ++itOut )
      {
      if( itIn.Get() != 0  && itIn.Get() != itCompare.Get() )
        {
        itOut.Set( itOut.Get() + 1 );
        std::cout << " Value at " << itOut.GetIndex()
                  << " set to " << itOut.Get()
                  << std::endl;
        }
      }
    }
  catch( itk::ExceptionObject& ex )
    {
    std::cout << "Exception:: " << __LINE__ << __FILE__ << std::endl;
    std::cout << ex.GetFile() << std::endl;
    std::cout << ex.GetLine() << std::endl;
    std::cout << ex.GetDescription() << std::endl;
    exit( EXIT_FAILURE );
    }
  return ambiguousMap;
}

LabelImagePointerType
BRAINSCutApplyModel
::CombineLabel( LabelImagePointerType& resultLabel,
                LabelImagePointerType& currentLabel,
                const unsigned char binaryToLabelValue )
{
  try
    {
    typedef itk::ImageRegionIterator<LabelImageType> RegionIteratorType;

    RegionIteratorType itOut( resultLabel, resultLabel->GetLargestPossibleRegion() );
    RegionIteratorType itIn( currentLabel, currentLabel->GetLargestPossibleRegion() );

    itIn.GoToBegin();
    itOut.GoToBegin();
    for( ; !itIn.IsAtEnd(); ++itIn, ++itOut )
      {
      if( !itIn.Get() )
        {
        continue;
        }
      if( binaryToLabelValue == 0 )
        {
        itOut.Set(itIn.Get() );
        }
      else
        {
        itOut.Set( binaryToLabelValue );
        }
      }
    }
  catch( itk::ExceptionObject& ex )
    {
    std::cout << "Exception:: " << __LINE__ << __FILE__ << std::endl;
    std::cout << ex.GetFile() << std::endl;
    std::cout << ex.GetLine() << std::endl;
    std::cout << ex.GetDescription() << std::endl;
    exit( EXIT_FAILURE );
    }
  return resultLabel;
}

LabelImagePointerType
BRAINSCutApplyModel
::PostProcessingANN( const std::string & continuousFilename,
                     scalarType threshold )
{
  WorkingImagePointer continuousImage = ReadImageByFilename( continuousFilename );

  LabelImagePointerType maskVolume;

  try
    {
    /* threshold */
    maskVolume = ThresholdImageAtLower( continuousImage, threshold);

    /* Get One label */
    maskVolume = Opening( maskVolume, 1 );
    maskVolume = GetOneConnectedRegion( maskVolume );

    /* opening and closing to get rid of island and holes */
    maskVolume = Closing( maskVolume, 1 );
    }
  catch( itk::ExceptionObject & ex )
    {
    std::cout << "ERROR at " << __LINE__ << "::" << __FILE__ << std::endl;
    std::cout << ex.what() << std::endl;
    }
  return maskVolume;
}

LabelImagePointerType
BRAINSCutApplyModel
::PostProcessingRF( const std::string & labelImageFilename,
                    const unsigned char labelValue )
{
  LabelImagePointerType labelImage =
    itkUtil::ReadImage<LabelImageType>( labelImageFilename );

  LabelImagePointerType resultBinaryImage;

  resultBinaryImage = FillHole(  labelImage, labelValue  );
  resultBinaryImage = Closing( resultBinaryImage, labelValue  );
  resultBinaryImage = GetOneConnectedRegion( resultBinaryImage );

  return resultBinaryImage;
}

LabelImagePointerType
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

  LabelImagePointerType mask = itkUtil::ScaleAndCast<WorkingImageType,
                                                     LabelImageType>(thresholder->GetOutput(),
                                                                     0,
                                                                     255);
  return mask;
}

WorkingImagePointer
BRAINSCutApplyModel
::ClipImageWithBinaryMask( WorkingImagePointer& image, WorkingImagePointer mask)
{
  std::cout << "Clip image ..." << std::endl;

  typedef itk::MultiplyImageFilter<WorkingImageType, WorkingImageType> ClipImageFilterType;
  ClipImageFilterType::Pointer clipper = ClipImageFilterType::New();
  clipper->SetInput1( image );
  clipper->SetInput2( mask );
  clipper->Update();

  return clipper->GetOutput();
}

LabelImagePointerType
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

  LabelImagePointerType mask = itkUtil::ScaleAndCast<WorkingImageType,
                                                     LabelImageType>(thresholder->GetOutput(),
                                                                     0,
                                                                     255);
  return mask;
}

LabelImagePointerType
BRAINSCutApplyModel
::ExtractLabel( const LabelImagePointerType& image, unsigned char thresholdValue  )
{
  typedef itk::BinaryThresholdImageFilter<LabelImageType, LabelImageType> ThresholdFilterType;
  ThresholdFilterType::Pointer thresholder = ThresholdFilterType::New();

  thresholder->SetInput( image );
  thresholder->SetInsideValue( 1 );
  thresholder->SetOutsideValue( 0 );
  thresholder->SetUpperThreshold( thresholdValue );
  thresholder->SetLowerThreshold( thresholdValue );
  thresholder->Update();

  LabelImagePointerType mask = itkUtil::TypeCast<LabelImageType, LabelImageType>( thresholder->GetOutput() );
  return mask;
}

LabelImagePointerType
BRAINSCutApplyModel
::Opening( LabelImagePointerType& image,
           const unsigned char labelValue )
{
  /*  Opening */
  // #include <itkBinaryOpeningByReconstructionImageFilter.h> consider this
  typedef itk::BinaryBallStructuringElement<LabelImageType::PixelType, DIMENSION> KernelType;
  KernelType           ball;
  KernelType::SizeType ballSize;

  /* Create the structuring element- a disk of radius 2 */
  ballSize.Fill(1);
  ball.SetRadius( ballSize );
  ball.CreateStructuringElement();
  typedef itk::BinaryMorphologicalOpeningImageFilter<LabelImageType,
                                                     LabelImageType,
                                                     KernelType> OpeningFilterType;
  OpeningFilterType::Pointer openingFilter = OpeningFilterType::New();
  openingFilter->SetInput( image );
  openingFilter->SetKernel( ball );
  openingFilter->SetForegroundValue( labelValue );
  openingFilter->Update();

  return openingFilter->GetOutput();
}

LabelImagePointerType
BRAINSCutApplyModel
::GetOneConnectedRegion( LabelImagePointerType& image)
{
  LabelImagePointerType resultMask;

  try
    {
    /* relabel images if they are disconnected */
    typedef itk::ConnectedComponentImageFilter<LabelImageType, LabelImageType>
      ConnectedBinaryImageFilterType;
    ConnectedBinaryImageFilterType::Pointer relabler = ConnectedBinaryImageFilterType::New();

    relabler->SetInput( image );

    /* relable images from the largest to smallest size */
    typedef itk::RelabelComponentImageFilter<LabelImageType, LabelImageType> RelabelInOrderFilterType;
    RelabelInOrderFilterType::Pointer relabelInOrder = RelabelInOrderFilterType::New();

    relabelInOrder->SetInput( relabler->GetOutput() );
    relabelInOrder->Update();

    LabelImagePointerType multipleLabelVolume = relabelInOrder->GetOutput();
    /* get the label one */
    constexpr unsigned char LargestLabelIndex = 1; //NOTE:  After RelabelInOrder 1 is always the largest remaining region.
    resultMask = ExtractLabel( multipleLabelVolume, LargestLabelIndex );
    }
  catch( itk::ExceptionObject& ex )
    {
    std::cout << "Exception:: " << __LINE__ << __FILE__ << std::endl;
    std::cout << ex.GetFile() << std::endl;
    std::cout << ex.GetLine() << std::endl;
    std::cout << ex.GetDescription() << std::endl;
    exit( EXIT_FAILURE );
    }

  return resultMask;
}

LabelImagePointerType
BRAINSCutApplyModel
::FillHole( LabelImagePointerType& mask,
            const unsigned char labelValue)
{
  typedef itk::BinaryFillholeImageFilter<LabelImageType> FillHoleFilterType;

  FillHoleFilterType::Pointer filler = FillHoleFilterType::New();
  filler->SetInput( mask );
  filler->SetFullyConnected( false );
  filler->SetForegroundValue( labelValue );
  filler->Update();
  return filler->GetOutput();
}

LabelImagePointerType
BRAINSCutApplyModel
::Closing( LabelImagePointerType& mask,
           const unsigned char labelValue )
{
  // NOTE: Consider this filter :#include <itkBinaryFillholeImageFilter.h>
  typedef itk::BinaryBallStructuringElement<LabelImageType::PixelType, DIMENSION> KernelType;
  KernelType           ball;
  KernelType::SizeType ballSize;

  /* Create the structuring element- a disk of radius 2 */
  std::cout << "Ball Size for Closing: 2 " << std::endl;
  ballSize.Fill(2);
  ball.SetRadius( ballSize );
  ball.CreateStructuringElement();

  /* Closing */
  typedef itk::BinaryMorphologicalClosingImageFilter<LabelImageType,
                                                     LabelImageType,
                                                     KernelType> ClosingFilterType;
  ClosingFilterType::Pointer closingFilter = ClosingFilterType::New();
  closingFilter->SetInput( mask );
  closingFilter->SetKernel( ball );
  closingFilter->SetForegroundValue( labelValue );
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

    cv::Mat openCVInputFeature = cv::Mat(1, inputVectorSize, CV_32F,
             &(it->second[0]) ); //
                                         // http://stackoverflow.com/questions/6701816/how-to-use-a-stdvector-in-a-c-function
    // std::cout<<FeatureInputVector::HashIndexFromKey( it->first );

    /* predict */
    if( m_method == "ANN" )
      {
      cv::Mat openCVOutput = cv::Mat::zeros(1, NumberOfROIS, CV_32S);
      std::cout<<openCVInputFeature<<std::endl;
      openCVOutput = this->m_openCVANN->predict( openCVInputFeature );

      /* insert result to the result output vector */
      resultOutputVector.insert( std::pair<hashKeyType, scalarType>(
                                 ( it->first ),
                                 openCVOutput.at<float>(0, roiNumber) ));
                                 //cvmGet( openCVOutput,  0, roiNumber) ) );
      // CV_MAT_ELEM( *openCVOutput, scalarType, 0, roiNumber) ) );
      }
    else if( m_method == "RandomForest" )
      {
      std::cout<<openCVInputFeature<<std::endl;
      const scalarType response = this->m_openCVRandomForest->predict( openCVInputFeature );
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
  m_computeSSE = sse;
  if( m_computeSSE )
    {
    m_ANNTestingSSEFileStream.open( this->m_myDataHandler->GetANNTestingSSEFilename().c_str(),
                                    std::fstream::out );
    if( !m_ANNTestingSSEFileStream.good() )
      {
      std::string errorMsg = " Cannot open the file! :";
      errorMsg += this->m_myDataHandler->GetANNTestingSSEFilename();
      std::cout << errorMsg << std::endl;
      throw BRAINSCutExceptionStringHandler( errorMsg );
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
  this->m_openCVANN = OpenCVMLPType::load( ANNModelFilename );
  //this->m_openCVANN->load( ANNModelFilename.c_str() );
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
  m_openCVRandomForest = cv::Algorithm::load<cv::ml::RTrees>( randomForestFilename.c_str() );
}

inline void
BRAINSCutApplyModel
::WritePredictROIProbabilityBasedOnReferenceImage( const PredictValueMapType& predictedOutput,
                                                   const WorkingImagePointer& referenceImage,
                                                   const WorkingImagePointer& roi,
                                                   const std::string & imageFilename,
                                                   const WorkingPixelType & labelValue )
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
::GetSubjectOutputDirectory( const DataSet& subject)
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
::GetContinuousPredictionFilename( const DataSet& subject, const std::string & currentROIName)
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
::GetROIVolumeName( const DataSet& subject, const std::string & currentROIName)
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

inline std::string
BRAINSCutApplyModel
::GetLabelMapFilename( const DataSet& subject )
{
  const std::string subjectID(subject.GetAttribute<StringValue>("Name") );

  const std::string & outputDir =  GetSubjectOutputDirectory( subject );
  const std::string & returnFilename = outputDir + "/" + subjectID + "_ANNLabel_seg.nii.gz";

  return returnFilename;
}

void
BRAINSCutApplyModel
::SetNumberOfTrees( const int trees)
{
  if( trees < 0 )
    {
    m_numberOfTrees = this->m_myDataHandler->GetMaxTreeCount();
    }
  else
    {
    m_numberOfTrees = trees;
    }
}

void
BRAINSCutApplyModel
::SetDepthOfTree( const int depth )
{
  if( depth < 0 )
    {
    m_depthOfTree = this->m_myDataHandler->GetMaxDepth();
    }
  else
    {
    m_depthOfTree = depth;
    }
}
