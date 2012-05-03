#ifndef FeatureInputVector_h
#define FeatureInputVector_h

#include "BRAINSCutDataHandler.h"

#include <itkLinearInterpolateImageFunction.h>
#include <itkGradientImageFilter.h>
/*
 * Note
 * - this class is to compute input vector of ONE subject for a given ROI
 *
 * in this class NO WARPING is occuring.
 * All the images should been warped onto subject space before hand
 */
class FeatureInputVector
{
public:
  FeatureInputVector();

  /** constants definition */
  static const scalarType MIN;
  static const scalarType MAX;

  /** type definition */
  typedef itk::GradientImageFilter<WorkingImageType,
                                   WorkingPixelType,
                                   WorkingPixelType> GradientFilterType;

  typedef itk::Image<itk::CovariantVector<WorkingPixelType, DIMENSION>, DIMENSION>::Pointer
    GradientImageType;

  typedef itk::Image<unsigned char, DIMENSION> BinaryImageType;

  typedef itk::LinearInterpolateImageFunction<WorkingImageType,
                                              WorkingPixelType> ImageLinearInterpolatorType;

  /* min/max */
  typedef std::pair<scalarType, scalarType> minmaxPairType;
  typedef std::vector<minmaxPairType>       minmaxPairVectorType;

  matrixType  GetInputVectorAt(WorkingImageVectorType & currentIndex);

  /** set functions */
  void SetGradientSize( unsigned int length);

  void SetImagesOfInterestInOrder( WorkingImageVectorType& images);

  void SetImagesOfSpatialLocation(  std::map<std::string, WorkingImagePointer>& SpatialLocationImages);

  void SetCandidateROIs( std::map<std::string, WorkingImagePointer>& candidateROIMap);

  void SetROIInOrder( DataSet::StringVectorType roiInOrder );

  void SetFeatureInputOfROI( std::map<std::string, WorkingImagePointer>& featureImages );

  void SetInputVectorSize();

  unsigned int GetInputVectorSize();

  void SetNormalization( const bool doNormalize);

  void NormalizationOfVector( InputVectorMapType& currentFeatureVector, std::string ROIName );

  /** get function(s) */
  InputVectorMapType GetFeatureInputOfROI( std::string ROIName );

  /* HashGenerator From Index */
  /* hash function is based on fixed size of 'size'
   * This would not work if the size of image bigger than the size we are using here.
   * The size could be easily changed though.
   */
  static int                         HashKeyFromIndex(const WorkingImageType::IndexType index);

  static WorkingImageType::IndexType HashIndexFromKey(const int offSet);

private:
  int          gradientSize;
  unsigned int inputVectorSize;
  bool         normalization;

  ImageLinearInterpolatorType::Pointer imageInterpolator;

  WorkingImageVectorType    imagesOfInterestInOrder;
  DataSet::StringVectorType roiIDsInOrder;

  /** deformed rho/phi/theta images*/
  std::map<std::string, WorkingImagePointer> spatialLocations;

  /** deformed candiateROIs */
  std::map<std::string, WorkingImagePointer> candidateROIs;

  /** gradient image of ROI */
  std::map<std::string, GradientImageType> gradientOfROI;

  /** feature output*/
  std::map<std::string, InputVectorMapType> featureInputOfROI;

  /** normalization parameters*/
  /*  mapping from ROIname to the vector of mean/max
   *  for given serios of imagesOfInterestInOrder
   */
  std::map<std::string, minmaxPairVectorType> minmax;

  /** private functions */
  void ComputeFeatureInputOfROI( std::string ROIName);

  void SetGradientImage( std::string ROIName );

  void SetNormalizationParameters( std::string ROIName);

  /** inline functions */
  inline void AddValueToElement( scalarType value, std::vector<scalarType>::iterator & elementIterator);

  inline void AddCandidateROIFeature( WorkingImageType::IndexType currentPixelIndex,
                                      std::vector<scalarType>::iterator & elementIterator);

  inline void AddSpatialLocation( WorkingImageType::IndexType currentPixelIndex,
                                  std::vector<scalarType>::iterator & elementIterator);

  inline void AddFeaturesImagesOfInterest( std::string ROIName, WorkingImageType::IndexType currentPixelIndex,
                                           std::vector<scalarType>::iterator & elementIterator);

  inline void AddFeaturesAlongGradient( std::string ROIName, WorkingImagePointer featureImage,
                                        WorkingImageType::IndexType currentPixelIndex,
                                        std::vector<scalarType>::iterator & elementIterator );

  inline std::map<std::string, scalarType> CalculateUnitDeltaAlongTheGradient(std::string ROIName,
                                                                              WorkingImageType::IndexType currentPixelIndex );

  inline std::pair<scalarType, scalarType>  SetMinMaxOfSubject( BinaryImageType::Pointer & labelImage,
                                                                const WorkingImagePointer & Image );
};
#endif
