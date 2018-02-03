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
#include "FeatureInputVector.h"
#include "BRAINSCutExceptionStringHandler.h"
#include "itkLabelStatisticsImageFilter.h"

constexpr unsigned int MAX_IMAGE_SIZE = 1024;
const WorkingImageType::IndexType ConstantHashIndexSize = {{1024, 1024, 1024}};

bool
FeatureInputVector
::DoScalingUnitTests()
{
  bool        allOK = true;
  const float x = 1.0F;
  const float min = 0.5F;
  const float max = 1.5F;
  // linear
  const float linearAnswer = 0.5F;

  if( linearAnswer != LinearScaling( x, min, max ) )
    {
    allOK = false;
    }

  if( !(doubleSigmoid( 0.5, 1, 1, 2) == Sigmoid( 0.5, 1, 1) ) )
    {
    allOK = false;
    }

  if( !(doubleSigmoid( 1.5, 1, 1, 2) == Sigmoid( 1.5, 1, 2) ) )
    {
    allOK = false;
    }

  if( !( ZScore( 0.0, 0, 1) != 0.0F ) )
    {
    allOK = false;
    }

  return allOK;
}

int
FeatureInputVector
::DoUnitTests(void) const
{
  std::cout << "INTERNAL TEST OF INDEX_KEY INPUTVECTOR" << std::endl;
  int allOK = EXIT_SUCCESS;

  for( WorkingImageType::IndexType::IndexValueType i = 0; i < ConstantHashIndexSize[0]; ++i )
    {
    for( WorkingImageType::IndexType::IndexValueType j = 0; j < ConstantHashIndexSize[1]; ++j )
      {
      for( WorkingImageType::IndexType::IndexValueType k = 0; k < ConstantHashIndexSize[2]; ++k )
        {
        // QuickTest
        WorkingImageType::IndexType index;
        index[0] = i;
        index[1] = j;
        index[2] = k;
        const hashKeyType                 currKey = HashKeyFromIndex( index );
        const WorkingImageType::IndexType outIndex = HashIndexFromKey(currKey);
        if( index != outIndex )
          {
          std::cout << "HACK: UNIT TEST " << index << " = " << outIndex << " with key " << currKey << std::endl;
          allOK = EXIT_FAILURE;
          }
        }
      }
    }
  if( allOK == EXIT_FAILURE )
    {
    std::cout << "ERROR: Hash lookups are not invertable." << std::endl;
    }
  else
    {
    std::cout << "All hash lookups for images < " << ConstantHashIndexSize << " are validated." << std::endl;
    }
  return allOK;
}

FeatureInputVector
::FeatureInputVector() :
  m_gradientSize(-1),
  m_inputVectorSize(0),
  m_normalizationMethod("None"),
  m_imageInterpolator(nullptr),
  m_imagesOfInterestInOrder(),
  m_roiIDsInOrder(),
  m_spatialLocations(),
  m_candidateROIs(),
  m_gradientOfROI(),
  m_minmax(),
  m_statistics(),
  m_featureNormalizationMap()
{
  m_spatialLocations.clear();
  m_candidateROIs.clear();
  m_gradientOfROI.clear();
  m_imageInterpolator = ImageLinearInterpolatorType::New();

  m_featureNormalizationMap["None"] = None;
  m_featureNormalizationMap["Linear"] = Linear;
  m_featureNormalizationMap["Sigmoid_Q05"] = Sigmoid_Q05;
  m_featureNormalizationMap["DoubleSigmoid_Q05"] = DoubleSigmoid_Q05;
  m_featureNormalizationMap["Sigmoid_Q01"] = Sigmoid_Q01;
  m_featureNormalizationMap["DoubleSigmoid_Q01"] = DoubleSigmoid_Q01;
  m_featureNormalizationMap["zScore"] = zScore;
  m_featureNormalizationMap["IQR"] = IQR;
}

FeatureInputVector
::~FeatureInputVector()
{
  this->m_imagesOfInterestInOrder.clear();
  this->m_spatialLocations.clear();
  this->m_gradientOfROI.clear();
  this->m_minmax.clear();
  this->m_candidateROIs.clear();
}

void
FeatureInputVector
::SetGradientSize(unsigned int length)
{
  m_gradientSize = length;
}

void
FeatureInputVector
::SetImagesOfInterestInOrder( WorkingImageVectorType& images)
{
  m_imagesOfInterestInOrder = images;
}

void
FeatureInputVector
::SetImagesOfSpatialLocation( std::map<std::string, WorkingImagePointer>& SpatialLocationImages)
{
  if( SpatialLocationImages.size() != 3 ||
      SpatialLocationImages.find("rho") == SpatialLocationImages.end() ||
      SpatialLocationImages.find("phi") == SpatialLocationImages.end()  ||
      SpatialLocationImages.find("theta") == SpatialLocationImages.end()  )
    {
    itkGenericExceptionMacro(<< "::number of images for spatial location should be 3 not "
                             << SpatialLocationImages.size() );
    }
  // Ensure that images are sufficiently small to process correctly.
  WorkingImageType::SizeType testSize =
    SpatialLocationImages.find("rho")->second->GetLargestPossibleRegion().GetSize();
  for( unsigned int q = 0; q < 3; ++q )
    {
    if( testSize[q] > MAX_IMAGE_SIZE )
      {
      std::cout << "ERROR: Image too large to process correctly" << std::endl;
      std::cout << testSize << " must be less than " << ConstantHashIndexSize << std::endl;
      exit(-1);
      }
    }
  m_spatialLocations = SpatialLocationImages;
}

void
FeatureInputVector
::SetCandidateROIs( std::map<std::string, WorkingImagePointer>& candidateROIMap)
{
  m_candidateROIs = candidateROIMap;
}

void
FeatureInputVector
::SetROIInOrder( DataSet::StringVectorType roiInOrder)
{
  m_roiIDsInOrder = roiInOrder;
}

void
FeatureInputVector
::SetGradientImage( std::string ROIName )
{
  GradientFilterType::Pointer gradientFilter = GradientFilterType::New();

  gradientFilter->SetInput( m_candidateROIs.find( ROIName)->second );
  try
    {
    gradientFilter->Update();
    }
  catch( ... )
    {
    std::string errorMsg = " Fail to generate itk gradient image.";
    throw BRAINSCutExceptionStringHandler( errorMsg );
    }
  m_gradientOfROI.insert( std::pair<std::string, GradientImageType>( ROIName, gradientFilter->GetOutput() ) );
}

void
FeatureInputVector
::SetInputVectorSize()
{
  if( m_candidateROIs.empty() || m_imagesOfInterestInOrder.empty() || m_spatialLocations.empty() )
    {
    std::string errorMsg = " Cannot compute input vector size properly.";
    errorMsg += "Either ROI(probability maps) or feature images has to be set to compute input vector size.";
    throw BRAINSCutExceptionStringHandler( errorMsg );
    }
  m_inputVectorSize = m_candidateROIs.size() + m_imagesOfInterestInOrder.size() * 3 + m_spatialLocations.size();
}

unsigned int
FeatureInputVector
::GetInputVectorSize()
{
  return m_inputVectorSize;
}

void
FeatureInputVector
::SetLinearNormalizationOn()
{
  SetNormalizationMethod( "Linear" );
}

void
FeatureInputVector
::SetNormalizationMethod( const std::string & normalizationMethod )
{
  switch( m_featureNormalizationMap[normalizationMethod] )
    {
    case None:
    case Linear:
    case Sigmoid_Q05:
    case DoubleSigmoid_Q05:
    case Sigmoid_Q01:
    case DoubleSigmoid_Q01:
    case zScore:
    case IQR:
      m_normalizationMethod = normalizationMethod;
      break;
    default:
      std::cout << "ERROR: No normalization method type support of "
                << normalizationMethod << std::endl;
      std::exit( EXIT_FAILURE );
    }
  return;
};
/*
InputVectorMapType
FeatureInputVector
::ComputeAndGetFeatureInputOfROI( std::string ROIName )
{
  std::map<std::string, InputVectorMapType> featureInputOfROI;
  if( featureInputOfROI.find( ROIName ) == featureInputOfROI.end() )
    {
    ComputeFeatureInputOfROI( ROIName);
    }
  return InputVectorMapType(featureInputOfROI.find( ROIName )->second);
}
*/
InputVectorMapType
FeatureInputVector
::ComputeAndGetFeatureInputOfROI( std::string ROIName)
{
  std::map<std::string, InputVectorMapType> featureInputOfROI;

  std::cout << "****************************************************" << std::endl;
  std::cout << "******** Compute Feature Input Of ROI **************" << std::endl;
  std::cout << "****************************************************" << std::endl;

  SetGradientImage( ROIName );

  typedef itk::ImageRegionIterator<WorkingImageType> ImageRegionIteratorType;

  InputVectorMapType currentFeatureVector;

  WorkingImagePointer currentROIImage = m_candidateROIs.find( ROIName)->second;

  /* iterate through each voxel in the probability map */
  ImageRegionIteratorType eachVoxelInROI( currentROIImage, currentROIImage->GetLargestPossibleRegion() );

  eachVoxelInROI.GoToBegin();

  while( !eachVoxelInROI.IsAtEnd() )
    {
    if( (eachVoxelInROI.Value() > (0.0F + FLOAT_TOLERANCE) ) &&
        (eachVoxelInROI.Value() < (1.0F - FLOAT_TOLERANCE) ) )
      {
      InputVectorType           oneRowInputFeature( m_inputVectorSize );
      InputVectorType::iterator featureElementIterator = oneRowInputFeature.begin();

      AddCandidateROIFeature( eachVoxelInROI.GetIndex(), featureElementIterator);
      AddSpatialLocation( eachVoxelInROI.GetIndex(), featureElementIterator);
      AddFeaturesImagesOfInterest(ROIName, eachVoxelInROI.GetIndex(), featureElementIterator);

      const int oneRowKey = FeatureInputVector::HashKeyFromIndex( eachVoxelInROI.GetIndex() );

      currentFeatureVector.insert( std::pair<hashKeyType, InputVectorType>( oneRowKey, oneRowInputFeature) );
      }
    ++eachVoxelInROI;
    }

  /* m_normalization */
  SetNormalizationParameters( ROIName );
  NormalizationOfVector( currentFeatureVector, ROIName );

  /* insert computed vector */
  featureInputOfROI.insert(std::pair<std::string, InputVectorMapType>( ROIName, currentFeatureVector) );

  return InputVectorMapType(featureInputOfROI.find( ROIName )->second);
}

/* set m_normalization parameters */
void
FeatureInputVector
::SetNormalizationParameters( std::string ROIName )
{
  constexpr unsigned char defaultLabel = 1;

  /* threshold roi */
  typedef itk::BinaryThresholdImageFilter<WorkingImageType,
                                          BinaryImageType> ThresholdType;

  ThresholdType::Pointer thresholder = ThresholdType::New();

  thresholder->SetInput( m_candidateROIs.find( ROIName)->second );
  thresholder->SetLowerThreshold( 0.0F + FLOAT_TOLERANCE );
  thresholder->SetInsideValue( defaultLabel );
  thresholder->Update();

  /* get min and max for each image type*/

  minmaxPairVectorType currentMinMaxVector;
  ImageTypeNo          ImageTypeNumber = 0;
  for( WorkingImageVectorType::const_iterator eachTypeOfImage = m_imagesOfInterestInOrder.begin();
       eachTypeOfImage != m_imagesOfInterestInOrder.end();
       ++eachTypeOfImage )
    {
    BinaryImageType::Pointer binaryImage = thresholder->GetOutput();
    minmaxPairType           eachMinMax = SetMinMaxOfSubject( binaryImage, *eachTypeOfImage );
    currentMinMaxVector.push_back( eachMinMax );

    // better do here
    typedef itk::LabelStatisticsImageFilter<WorkingImageType, BinaryImageType> StatisticCalculatorType;
    StatisticCalculatorType::Pointer statisticCalculator = StatisticCalculatorType::New();

    statisticCalculator->SetInput( *eachTypeOfImage );
    statisticCalculator->SetLabelInput( binaryImage );

    statisticCalculator->SetHistogramParameters( 100,
                                                 0.0F,
                                                 1.0F );
    statisticCalculator->SetUseHistograms( true);
    statisticCalculator->Update();
    m_statistics[ROIName][ImageTypeNumber]["Minimum"] =  statisticCalculator->GetMinimum( defaultLabel );
    m_statistics[ROIName][ImageTypeNumber]["Maximum"] =  statisticCalculator->GetMaximum( defaultLabel );
    m_statistics[ROIName][ImageTypeNumber]["Mean"] =  statisticCalculator->GetMean( defaultLabel );
    m_statistics[ROIName][ImageTypeNumber]["Sigma"] =  statisticCalculator->GetSigma( defaultLabel );
    m_statistics[ROIName][ImageTypeNumber]["Median"] =  statisticCalculator->GetMedian( defaultLabel );

    StatisticCalculatorType::HistogramPointer histogram = statisticCalculator->GetHistogram( defaultLabel );
    m_statistics[ROIName][ImageTypeNumber]["Q_01"] = histogram->Quantile( 0, 0.01 );
    m_statistics[ROIName][ImageTypeNumber]["Q_99"] = histogram->Quantile( 0, 0.99 );
    m_statistics[ROIName][ImageTypeNumber]["Q_25"] = histogram->Quantile( 0, 0.25 );
    m_statistics[ROIName][ImageTypeNumber]["Q_75"] = histogram->Quantile( 0, 0.75 );
    m_statistics[ROIName][ImageTypeNumber]["Q_95"] = histogram->Quantile( 0, 0.95 );
    m_statistics[ROIName][ImageTypeNumber]["Q_05"] = histogram->Quantile( 0, 0.05 );

    std::cout << ROIName << " :: " << ImageTypeNumber << " :: Minimum :: "
              << m_statistics[ROIName][ImageTypeNumber]["Minimum"]  << std::endl;
    std::cout << ROIName << " :: " << ImageTypeNumber << " :: Maximum :: "
              << m_statistics[ROIName][ImageTypeNumber]["Maximum"]  << std::endl;
    std::cout << ROIName << " :: " << ImageTypeNumber << " :: Mean :: "
              << m_statistics[ROIName][ImageTypeNumber]["Mean"]  << std::endl;
    std::cout << ROIName << " :: " << ImageTypeNumber << " :: Sigma :: "
              << m_statistics[ROIName][ImageTypeNumber]["Sigma"]  << std::endl;
    std::cout << ROIName << " :: " << ImageTypeNumber << " :: Median :: "
              << m_statistics[ROIName][ImageTypeNumber]["Median"]  << std::endl;
    std::cout << ROIName << " :: " << ImageTypeNumber << " :: Q_01 :: "
              << m_statistics[ROIName][ImageTypeNumber]["Q_01"]  << std::endl;
    std::cout << ROIName << " :: " << ImageTypeNumber << " :: Q_05 :: "
              << m_statistics[ROIName][ImageTypeNumber]["Q_05"]  << std::endl;
    std::cout << ROIName << " :: " << ImageTypeNumber << " :: Q_25 :: "
              << m_statistics[ROIName][ImageTypeNumber]["Q_25"]  << std::endl;
    std::cout << ROIName << " :: " << ImageTypeNumber << " :: Q_75 :: "
              << m_statistics[ROIName][ImageTypeNumber]["Q_75"]  << std::endl;
    std::cout << ROIName << " :: " << ImageTypeNumber << " :: Q_95 :: "
              << m_statistics[ROIName][ImageTypeNumber]["Q_95"]  << std::endl;
    std::cout << ROIName << " :: " << ImageTypeNumber << " :: Q_99 :: "
              << m_statistics[ROIName][ImageTypeNumber]["Q_99"]  << std::endl;
    ++ImageTypeNumber;
    }

  m_minmax[ROIName] = currentMinMaxVector;
}

/** inline functions */
inline void
FeatureInputVector
::AddValueToElement( scalarType value, std::vector<scalarType>::iterator & elementIterator)
{
  try
    {
    *elementIterator = value;
    ++elementIterator;
    }
  catch( ... )
    {
    std::string errorMsg = "Fail To Add Value To Element.";
    throw BRAINSCutExceptionStringHandler( errorMsg );
    }
  // std::cout<<value<<" ";
}

inline void
FeatureInputVector
::AddCandidateROIFeature( WorkingImageType::IndexType currentPixelIndex,
                          std::vector<scalarType>::iterator & elementIterator)
{
  for( DataSet::StringVectorType::const_iterator roiStringIt = m_roiIDsInOrder.begin();
       roiStringIt != m_roiIDsInOrder.end();
       ++roiStringIt )  // iterate each ROI candidates in order specified in "roi IDs in order"
    {
    WorkingPixelType currentProbability = m_candidateROIs.find( *roiStringIt )->second->GetPixel( currentPixelIndex );
    if( currentProbability > 0.0F +  FLOAT_TOLERANCE )
      {
      AddValueToElement( HundredPercentValue, elementIterator );
      }
    else
      {
      AddValueToElement( ZeroPercentValue, elementIterator );
      }
    }
}

inline void
FeatureInputVector
::AddSpatialLocation( WorkingImageType::IndexType currentPixelIndex,
                      std::vector<scalarType>::iterator & elementIterator)
{
  // std::cout<<" (spatial) ";
  AddValueToElement( m_spatialLocations.find("rho")->second->GetPixel( currentPixelIndex ), elementIterator );
  AddValueToElement( m_spatialLocations.find("phi")->second->GetPixel( currentPixelIndex ), elementIterator );
  AddValueToElement( m_spatialLocations.find("theta")->second->GetPixel( currentPixelIndex ), elementIterator );
}

inline void
FeatureInputVector
::AddFeaturesImagesOfInterest(  std::string ROIName,
                                WorkingImageType::IndexType currentPixelIndex,
                                std::vector<scalarType>::iterator & elementIterator )
{
  for( WorkingImageVectorType::const_iterator wit = m_imagesOfInterestInOrder.begin();
       wit != m_imagesOfInterestInOrder.end();
       ++wit )
    {
    // std::cout<<"(IMG) ";
    AddFeaturesAlongGradient( ROIName, (*wit), currentPixelIndex, elementIterator);
    }
}

inline void
FeatureInputVector
::AddFeaturesAlongGradient( std::string ROIName,
                            const WorkingImagePointer& featureImage,
                            WorkingImageType::IndexType currentPixelIndex,
                            std::vector<scalarType>::iterator & elementIterator )
{
  const std::map<std::string, scalarType> delta = CalculateUnitDeltaAlongTheGradient( ROIName,
                                                                                      currentPixelIndex );
  itk::Point<WorkingPixelType, DIMENSION> CenterPhysicalPoint;

  featureImage->TransformIndexToPhysicalPoint( currentPixelIndex, CenterPhysicalPoint );

  m_imageInterpolator->SetInputImage( featureImage );
  const float gradientUnitSize = 1.0F;
  for( float i = -m_gradientSize; i <= m_gradientSize; i = i + gradientUnitSize )
    {
    // std::cout<<"(gradient at "<< i<<" )";
    itk::Point<WorkingPixelType, 3> gradientLocation = CenterPhysicalPoint;

    gradientLocation[0] = CenterPhysicalPoint[0] + i * (delta.find("deltaX")->second);
    gradientLocation[1] = CenterPhysicalPoint[1] + i * (delta.find("deltaY")->second);
    gradientLocation[2] = CenterPhysicalPoint[2] + i * (delta.find("deltaZ")->second);

    itk::ContinuousIndex<WorkingPixelType, DIMENSION> ContinuousIndexOfGradientLocation;

    featureImage->TransformPhysicalPointToContinuousIndex( gradientLocation, ContinuousIndexOfGradientLocation );

    AddValueToElement( static_cast<scalarType>( m_imageInterpolator->
                                                EvaluateAtContinuousIndex( ContinuousIndexOfGradientLocation ) ),
                       elementIterator );
    }
}

inline std::map<std::string, scalarType>
FeatureInputVector
::CalculateUnitDeltaAlongTheGradient( std::string ROIName,
                                      WorkingImageType::IndexType currentPixelIndex )
{
  WorkingPixelType deltaX = m_gradientOfROI.find( ROIName )->second->GetPixel( currentPixelIndex )[0];
  WorkingPixelType deltaY = m_gradientOfROI.find( ROIName )->second->GetPixel( currentPixelIndex )[1];
  WorkingPixelType deltaZ = m_gradientOfROI.find( ROIName )->second->GetPixel( currentPixelIndex )[2];

  const scalarType Length = std::sqrt(deltaX * deltaX + deltaY * deltaY + deltaZ * deltaZ);
  const scalarType inverseLength =  ( Length > 0.0F ) ? 1.0 / Length : 1;

  std::map<std::string, scalarType> unitGradient;

  unitGradient["deltaX"] = deltaX * inverseLength;
  unitGradient["deltaY"] = deltaY * inverseLength;
  unitGradient["deltaZ"] = deltaZ * inverseLength;

  return unitGradient;
}

/* Hash Generator from index */
hashKeyType
FeatureInputVector
::HashKeyFromIndex( const WorkingImageType::IndexType index )
{
  /*
   * calculating offset
   * hashValue = i[2] + i[1]*s[1] + i[0]*s[0]*s[1]
   */
  hashKeyType hashValue = 0; // TODO HACK REGINA: HashKeys should be unsigned!

  const unsigned int lastDimensionIndex = DIMENSION - 1;

  for( unsigned int i = 0; i < (lastDimensionIndex); ++i )
    {
    hashValue += index[i];
    hashValue *= ConstantHashIndexSize[i];
    }
  hashValue += index[lastDimensionIndex];
  return hashValue;
}

WorkingImageType::IndexType
FeatureInputVector
::HashIndexFromKey(const hashKeyType offSet) // TODO HACK REGINA: HashKeys should be unsigned!
{
  WorkingImageType::IndexType key;
  hashKeyType                 remainedOffSet = offSet; // TODO HACK REGINA: HashKeys should be unsigned!

  for( int d = DIMENSION - 1; d >= 0; d-- )
    {
    key[d] = remainedOffSet % ConstantHashIndexSize[d];
    remainedOffSet = remainedOffSet / ConstantHashIndexSize[d];
    }
  return key;
}

inline std::pair<scalarType, scalarType>
FeatureInputVector
::SetMinMaxOfSubject( BinaryImageType::Pointer & labelImage, const WorkingImagePointer& Image )
{
  typedef itk::LabelStatisticsImageFilter<WorkingImageType, BinaryImageType> StatisticCalculatorType;
  StatisticCalculatorType::Pointer statisticCalculator = StatisticCalculatorType::New();

  statisticCalculator->SetInput( Image );
  statisticCalculator->SetLabelInput( labelImage);

  statisticCalculator->Update();

  /*std::cout << " * Min : " << statisticCalculator->GetMinimum(1)
            << " * Max : " << statisticCalculator->GetMaximum(1)
            << std::endl; */
  return minmaxPairType( std::pair<scalarType, scalarType>( statisticCalculator->GetMinimum(1),
                                                            statisticCalculator->GetMaximum(1) ) );
}

void
FeatureInputVector
::NormalizationOfVector( InputVectorMapType& currentFeatureVector, std::string ROIName )
{
  for( InputVectorMapType::iterator eachInputVector = currentFeatureVector.begin();
       eachInputVector != currentFeatureVector.end();
       ++eachInputVector )
    {
    InputVectorType::iterator featureElementIterator = (eachInputVector->second).begin();
    featureElementIterator += (m_roiIDsInOrder.size() + m_spatialLocations.size() );
    /*
    minmaxPairVectorType currentMinmaxPairVector = m_minmax.find(ROIName)->second;
      for( minmaxPairVectorType::const_iterator m_minmaxIt = currentMinmaxPairVector.begin();
         m_minmaxIt != currentMinmaxPairVector.end();
         ++m_minmaxIt )
      {
      for(  float i = -m_gradientSize; i <= m_gradientSize; i = i + 1.0F )
        {
        if( m_normalizationMethod == "Linear" )
            {
            scalarType normalizedValue = ( *featureElementIterator - m_minmaxIt->first );
            normalizedValue = normalizedValue / ( m_minmaxIt->second - m_minmaxIt->first );
            *featureElementIterator = normalizedValue;
            }
        ++featureElementIterator;
        }
      }
    */
    ImageTypeNo currentImgType = 0;
    for( WorkingImageVectorType::const_iterator eachTypeOfImage = m_imagesOfInterestInOrder.begin();
         eachTypeOfImage != m_imagesOfInterestInOrder.end();
         ++eachTypeOfImage )
      {
      for(  float i = -m_gradientSize; i <= m_gradientSize; i = i + 1.0F )
        {
        if( m_normalizationMethod == "Linear" )
          {
          *featureElementIterator = LinearScaling( *featureElementIterator,
                                                   m_statistics[ROIName][currentImgType]["Minimum"],
                                                   m_statistics[ROIName][currentImgType]["Maximum"] );
          }
        else if( m_normalizationMethod == "Sigmoid_Q05" )
          {
          // Sigmoid
          *featureElementIterator = Sigmoid( *featureElementIterator,
                                             m_statistics[ROIName][currentImgType]["Median"],
                                             m_statistics[ROIName][currentImgType]["Q_95"]
                                             - m_statistics[ROIName][currentImgType]["Q_05"] );
          }
        else if( m_normalizationMethod == "Sigmoid_Q01" )
          {
          // Sigmoid
          *featureElementIterator = Sigmoid( *featureElementIterator,
                                             m_statistics[ROIName][currentImgType]["Median"],
                                             m_statistics[ROIName][currentImgType]["Q_99"]
                                             - m_statistics[ROIName][currentImgType]["Q_01"] );
          }
        else if( m_normalizationMethod == "DoubleSigmoid_Q05" )
          {
          *featureElementIterator = doubleSigmoid( *featureElementIterator,
                                                   m_statistics[ROIName][currentImgType]["Median"],
                                                   m_statistics[ROIName][currentImgType]["Median"]
                                                   - m_statistics[ROIName][currentImgType]["Q_05"],
                                                   m_statistics[ROIName][currentImgType]["Q_95"]
                                                   - m_statistics[ROIName][currentImgType]["Median"] );
          }
        else if( m_normalizationMethod == "DoubleSigmoid_Q01" )
          {
          *featureElementIterator = doubleSigmoid( *featureElementIterator,
                                                   m_statistics[ROIName][currentImgType]["Median"],
                                                   m_statistics[ROIName][currentImgType]["Median"]
                                                   - m_statistics[ROIName][currentImgType]["Q_01"],
                                                   m_statistics[ROIName][currentImgType]["Q_99"]
                                                   - m_statistics[ROIName][currentImgType]["Median"] );
          }
        else if( m_normalizationMethod == "zScore" )
          {
          *featureElementIterator = ZScore( *featureElementIterator,
                                            m_statistics[ROIName][currentImgType]["Mean"],
                                            m_statistics[ROIName][currentImgType]["Sigma"] );
          }
        else if( m_normalizationMethod == "IQR" )
          {
          float current_IQR = m_statistics[ROIName][currentImgType]["Q_75"]
            - m_statistics[ROIName][currentImgType]["Q_25"];
          *featureElementIterator = ZScore( *featureElementIterator,
                                            m_statistics[ROIName][currentImgType]["Median"],
                                            current_IQR );
          }
        else if( m_normalizationMethod == "None" )
          {
          // do nothing
          }
        else
          {
          std::cout << "In valid normalization type of " << m_normalizationMethod << std::endl;
          std::exit( EXIT_FAILURE);
          }
        ++featureElementIterator;
        }
      ++currentImgType;
      }
    }
}

/* normalization is all to zero to one range */
float
FeatureInputVector
::Sigmoid( const float x,
           const float center,
           const float range)
{
  //    1
  // _____________
  // 1 + exp( -2x )
  float exp_value;
  float return_value;

  /*** Exponential calculation ***/
  const float fixedSlopeConstant = 8.0F;

  exp_value = exp( (double) -fixedSlopeConstant * ( x - center )  / range);

  /*** Final Sigmoid value ***/
  float sigmoid_value = 1 / (1 + exp_value);
  return_value =  LinearScaling( sigmoid_value, 0.0F, 1.0F);
  // std::cout<< " x = " << x << ", "
  //         << " center = " << center << ", "
  //         << " range = " << range << ", "
  //         << " sigmoid = " << sigmoid_value <<","
  //         << " returnValue = "<< return_value << std::endl;

  return return_value;
}

float
FeatureInputVector
::LinearScaling( const float x,
                 const float min,
                 const float max )
{
  if( min == ZeroPercentValue && max == HundredPercentValue )
    {
    return x;
    }

  float newRange = HundredPercentValue - ZeroPercentValue;
  float oldRange = max - min;
  float return_value = ZeroPercentValue + ( x - min ) * ( newRange / oldRange );

  return return_value;
}

float
FeatureInputVector
::ZScore( const float x,
          const float mu,
          const float sigma )
{
  float return_value;
  float zScore_value = ( x - mu ) / sigma;

  float sigma4Range = 4.0F;

  return_value = LinearScaling( zScore_value, -sigma4Range, sigma4Range );
  // std::cout<< " x = " << x << ", "
  //          << " mu = " << mu << ", "
  //          << " sigma = " << sigma << ", "
  //          << " zScore = " << zScore_value<<","
  //          << " returnValue = "<< return_value << std::endl;
  return return_value;
}

float
FeatureInputVector
::doubleSigmoid( const float x,
                 const float t,
                 const float r1,
                 const float r2)
{
  //           1
  // ________________________
  // 1 + exp( -2 * (s-t)/r )
  float exp_value;
  float return_value;
  float r;

  if( x > t )
    {
    r = r2;
    }
  else
    {
    r = r1;
    }

  /*** Exponential calculation ***/
  const float fixedSlopeConstant = 8.0F;
  exp_value = exp( (double) -fixedSlopeConstant * ( x - t ) / r);

  /*** Final Sigmoid value ***/
  float sigmoid_value = 1 / (1 + exp_value);

  return_value = LinearScaling( sigmoid_value, 0.0F, 1.0F);

  // std::cout<< " x = " << x << ", "
  //         << " t = " << t << ", "
  //         << " r1 = " << r1 << ", "
  //         << " r2 = " << r2 << ", "
  //         << " sigmoid = " << sigmoid_value <<","
  //         << " returnValue = "<< return_value << std::endl;

  return return_value;
}
