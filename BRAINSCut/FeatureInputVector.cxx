#include "FeatureInputVector.h"
#include "BRAINSCutExceptionStringHandler.h"
#include "itkLabelStatisticsImageFilter.h"

const scalarType FeatureInputVector::MIN = -1.0F;
const scalarType FeatureInputVector::MAX = 1.0F;

const unsigned int                MAX_IMAGE_SIZE = 1024;
const WorkingImageType::IndexType ConstantHashIndexSize = {{1024, 1024, 1024}};

int
FeatureInputVector
::DoUnitTests(void) const
{
  std::cout << "INTERNAL TEST OF INDEX_KEY INPUTVECTOR" << std::endl;
  int allOK = EXIT_SUCCESS;

  for( WorkingImageType::IndexType::IndexValueType i = 0; i < ConstantHashIndexSize[0]; i++ )
    {
    for( WorkingImageType::IndexType::IndexValueType j = 0; j < ConstantHashIndexSize[1]; j++ )
      {
      for( WorkingImageType::IndexType::IndexValueType k = 0; k < ConstantHashIndexSize[2]; k++ )
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
  gradientSize(-1),
  normalization(false)
{
  spatialLocations.clear();
  candidateROIs.clear();
  gradientOfROI.clear();
  featureInputOfROI.clear();
  imageInterpolator = ImageLinearInterpolatorType::New();
}

void
FeatureInputVector
::SetGradientSize(unsigned int length)
{
  gradientSize = length;
}

void
FeatureInputVector
::SetImagesOfInterestInOrder( WorkingImageVectorType& images)
{
  imagesOfInterestInOrder = images;
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
  for( unsigned int q = 0; q < 3; q++ )
    {
    if( testSize[q] > MAX_IMAGE_SIZE )
      {
      std::cout << "ERROR: Image too large to process correctly" << std::endl;
      std::cout << testSize << " must be less than " << ConstantHashIndexSize << std::endl;
      exit(-1);
      }
    }
  spatialLocations = SpatialLocationImages;
}

void
FeatureInputVector
::SetCandidateROIs( std::map<std::string, WorkingImagePointer>& candidateROIMap)
{
  candidateROIs = candidateROIMap;
}

void
FeatureInputVector
::SetROIInOrder( DataSet::StringVectorType roiInOrder)
{
  roiIDsInOrder = roiInOrder;
}

void
FeatureInputVector
::SetGradientImage( std::string ROIName )
{
  GradientFilterType::Pointer gradientFilter = GradientFilterType::New();

  gradientFilter->SetInput( candidateROIs.find( ROIName)->second );
  try
    {
    gradientFilter->Update();
    }
  catch( ... )
    {
    std::string errorMsg = " Fail to generate itk gradient image.";
    throw BRAINSCutExceptionStringHandler( errorMsg );
    }
  gradientOfROI.insert( std::pair<std::string, GradientImageType>( ROIName, gradientFilter->GetOutput() ) );
}

void
FeatureInputVector
::SetInputVectorSize()
{
  if( candidateROIs.empty() || imagesOfInterestInOrder.empty() || spatialLocations.empty() )
    {
    std::string errorMsg = " Cannot compute input vector size properly.";
    errorMsg += "Either ROI(probability maps) or feature images has to be set to compute input vector size.";
    throw BRAINSCutExceptionStringHandler( errorMsg );
    }
  inputVectorSize = candidateROIs.size() + imagesOfInterestInOrder.size() * 3 + spatialLocations.size();
}

unsigned int
FeatureInputVector
::GetInputVectorSize()
{
  return inputVectorSize;
}

void
FeatureInputVector
::SetNormalization( const bool doNormalize)
{
  normalization = doNormalize;
};

InputVectorMapType
FeatureInputVector
::GetFeatureInputOfROI( std::string ROIName )
{
  if( featureInputOfROI.find( ROIName ) == featureInputOfROI.end() )
    {
    ComputeFeatureInputOfROI( ROIName);
    }
  return InputVectorMapType(featureInputOfROI.find( ROIName )->second);
}

void
FeatureInputVector
::ComputeFeatureInputOfROI( std::string ROIName)
{
  SetGradientImage( ROIName );

  typedef itk::ImageRegionIterator<WorkingImageType> ImageRegionIteratorType;

  InputVectorMapType currentFeatureVector;

  WorkingImagePointer currentROIImage = candidateROIs.find( ROIName)->second;

  /* iterate through each voxel in the probability map */
  ImageRegionIteratorType eachVoxelInROI( currentROIImage, currentROIImage->GetLargestPossibleRegion() );

  eachVoxelInROI.GoToBegin();

  while( !eachVoxelInROI.IsAtEnd() )
    {
    if( (eachVoxelInROI.Value() > (0.0F + FLOAT_TOLERANCE) ) && (eachVoxelInROI.Value() < (1.0F - FLOAT_TOLERANCE) ) )
      {
      InputVectorType           oneRowInputFeature( inputVectorSize );
      InputVectorType::iterator featureElementIterator = oneRowInputFeature.begin();

      AddCandidateROIFeature( eachVoxelInROI.GetIndex(), featureElementIterator);
      AddSpatialLocation( eachVoxelInROI.GetIndex(), featureElementIterator);
      AddFeaturesImagesOfInterest(ROIName, eachVoxelInROI.GetIndex(), featureElementIterator);

      const int oneRowKey = FeatureInputVector::HashKeyFromIndex( eachVoxelInROI.GetIndex() );

      currentFeatureVector.insert( std::pair<hashKeyType, InputVectorType>( oneRowKey, oneRowInputFeature) );
      }
    ++eachVoxelInROI;
    }

  /* normalization */
  if( normalization )
    {
    SetNormalizationParameters( ROIName );
    NormalizationOfVector( currentFeatureVector, ROIName );
    }

  /* insert computed vector */
  featureInputOfROI.insert(std::pair<std::string, InputVectorMapType>( ROIName, currentFeatureVector) );
}

/* set normalization parameters */
void
FeatureInputVector
::SetNormalizationParameters( std::string ROIName )
{
  /* threshold roi */
  typedef itk::BinaryThresholdImageFilter<WorkingImageType,
                                          BinaryImageType> ThresholdType;

  ThresholdType::Pointer thresholder = ThresholdType::New();

  thresholder->SetInput( candidateROIs.find( ROIName)->second );
  thresholder->SetLowerThreshold( 0.0F + FLOAT_TOLERANCE );
  thresholder->SetInsideValue(1);
  thresholder->Update();

  /* get min and max for each image type*/

  minmaxPairVectorType currentMinMaxVector;
  for( WorkingImageVectorType::const_iterator eachTypeOfImage = imagesOfInterestInOrder.begin();
       eachTypeOfImage != imagesOfInterestInOrder.end();
       ++eachTypeOfImage )
    {
    BinaryImageType::Pointer binaryImage = thresholder->GetOutput();
    minmaxPairType           eachMinMax = SetMinMaxOfSubject( binaryImage, *eachTypeOfImage );
    currentMinMaxVector.push_back( eachMinMax );
    }

  minmax[ROIName] = currentMinMaxVector;
}

/** inline functions */
inline void
FeatureInputVector
::AddValueToElement( scalarType value, std::vector<scalarType>::iterator & elementIterator)
{
  try
    {
    *elementIterator = value;
    elementIterator++;
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
  for( DataSet::StringVectorType::const_iterator roiStringIt = roiIDsInOrder.begin();
       roiStringIt != roiIDsInOrder.end();
       ++roiStringIt )  // iterate each ROI candidates in order specified in "roi IDs in order"
    {
    WorkingPixelType currentProbability = candidateROIs.find( *roiStringIt )->second->GetPixel( currentPixelIndex );
    if( currentProbability > 0.0F +  FLOAT_TOLERANCE )
      {
      AddValueToElement( MAX, elementIterator );
      }
    else
      {
      AddValueToElement( MIN, elementIterator );
      }
    }
}

inline void
FeatureInputVector
::AddSpatialLocation( WorkingImageType::IndexType currentPixelIndex,
                      std::vector<scalarType>::iterator & elementIterator)
{
  // std::cout<<" (spatial) ";
  AddValueToElement( spatialLocations.find("rho")->second->GetPixel( currentPixelIndex ), elementIterator );
  AddValueToElement( spatialLocations.find("phi")->second->GetPixel( currentPixelIndex ), elementIterator );
  AddValueToElement( spatialLocations.find("theta")->second->GetPixel( currentPixelIndex ), elementIterator );
}

inline void
FeatureInputVector
::AddFeaturesImagesOfInterest(  std::string ROIName,
                                WorkingImageType::IndexType currentPixelIndex,
                                std::vector<scalarType>::iterator & elementIterator )
{
  for( WorkingImageVectorType::const_iterator wit = imagesOfInterestInOrder.begin();
       wit != imagesOfInterestInOrder.end();
       ++wit )
    {
    // std::cout<<"(IMG) ";
    AddFeaturesAlongGradient( ROIName, (*wit), currentPixelIndex, elementIterator);
    }
}

inline void
FeatureInputVector
::AddFeaturesAlongGradient( std::string ROIName,
                            WorkingImagePointer featureImage,
                            WorkingImageType::IndexType currentPixelIndex,
                            std::vector<scalarType>::iterator & elementIterator )
{
  const std::map<std::string, scalarType> delta = CalculateUnitDeltaAlongTheGradient( ROIName,
                                                                                      currentPixelIndex );
  itk::Point<WorkingPixelType, DIMENSION> CenterPhysicalPoint;

  featureImage->TransformIndexToPhysicalPoint( currentPixelIndex, CenterPhysicalPoint );

  imageInterpolator->SetInputImage( featureImage );
  for( float i = -gradientSize; i <= gradientSize; i = i + 1.0F )
    {
    // std::cout<<"(gradient at "<< i<<" )";
    itk::Point<WorkingPixelType, 3> gradientLocation = CenterPhysicalPoint;

    gradientLocation[0] = CenterPhysicalPoint[0] + i * (delta.find("deltaX")->second);
    gradientLocation[1] = CenterPhysicalPoint[1] + i * (delta.find("deltaY")->second);
    gradientLocation[2] = CenterPhysicalPoint[2] + i * (delta.find("deltaZ")->second);

    itk::ContinuousIndex<WorkingPixelType, DIMENSION> ContinuousIndexOfGradientLocation;

    featureImage->TransformPhysicalPointToContinuousIndex( gradientLocation, ContinuousIndexOfGradientLocation );

    AddValueToElement( static_cast<scalarType>( imageInterpolator->
                                                EvaluateAtContinuousIndex( ContinuousIndexOfGradientLocation ) ),
                       elementIterator );
    }
}

inline std::map<std::string, scalarType>
FeatureInputVector
::CalculateUnitDeltaAlongTheGradient( std::string ROIName,
                                      WorkingImageType::IndexType currentPixelIndex )
{
  WorkingPixelType deltaX = gradientOfROI.find( ROIName )->second->GetPixel( currentPixelIndex )[0];
  WorkingPixelType deltaY = gradientOfROI.find( ROIName )->second->GetPixel( currentPixelIndex )[1];
  WorkingPixelType deltaZ = gradientOfROI.find( ROIName )->second->GetPixel( currentPixelIndex )[2];

  const scalarType Length = vcl_sqrt(deltaX * deltaX + deltaY * deltaY + deltaZ * deltaZ);
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

  for( unsigned int i = 0; i < (lastDimensionIndex); i++ )
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
  minmaxPairVectorType currentMinmaxPairVector = minmax.find(ROIName)->second;

  for( InputVectorMapType::iterator eachInputVector = currentFeatureVector.begin();
       eachInputVector != currentFeatureVector.end();
       ++eachInputVector )
    {
    InputVectorType::iterator featureElementIterator = (eachInputVector->second).begin();
    featureElementIterator += (roiIDsInOrder.size() + spatialLocations.size() );
    for( minmaxPairVectorType::const_iterator minmaxIt = currentMinmaxPairVector.begin();
         minmaxIt != currentMinmaxPairVector.end();
         ++minmaxIt )
      {
      for(  float i = -gradientSize; i <= gradientSize; i = i + 1.0F )
        {
        scalarType normalizedValue = ( *featureElementIterator - minmaxIt->first );
        normalizedValue = normalizedValue / ( minmaxIt->second - minmaxIt->first );
        *featureElementIterator = normalizedValue;
        featureElementIterator++;
        }
      }
    }
}
