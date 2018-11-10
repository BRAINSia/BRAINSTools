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
#ifndef FeatureInputVector_h
#define FeatureInputVector_h

#include "BRAINSCutDataHandler.h"

#include <itkLinearInterpolateImageFunction.h>
#include <itkGradientImageFilter.h>

using hashKeyType = unsigned int;
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
  ~FeatureInputVector();

  int DoUnitTests() const; // A series of unit tests to verify results.

  bool DoScalingUnitTests();

  /* normalization type */
  enum FeatureNormalizationMethodEnum
    {
    None,
    Linear,
    Sigmoid_Q05,
    DoubleSigmoid_Q05,
    Sigmoid_Q01,
    DoubleSigmoid_Q01,
    zScore,
    IQR
    };

  /** type definition */
  using GradientFilterType = itk::GradientImageFilter<WorkingImageType,
                                   WorkingPixelType,
                                   WorkingPixelType>;

  typedef itk::Image<itk::CovariantVector<WorkingPixelType, DIMENSION>, DIMENSION>::Pointer
    GradientImageType;

  using BinaryImageType = itk::Image<unsigned char, DIMENSION>;

  using ImageLinearInterpolatorType = itk::LinearInterpolateImageFunction<WorkingImageType,
                                              WorkingPixelType>;

  /* min/max */
  using minmaxPairType = std::pair<scalarType, scalarType>;
  using minmaxPairVectorType = std::vector<minmaxPairType>;

  /* normalizationParameters */
  using ImageTypeNo = unsigned int;
  using ROITypeString = std::string;
  using StatisticsString = std::string;

  using NormalizationParameterType = std::map<ROITypeString,
                   std::map<ImageTypeNo, std::map<StatisticsString, scalarType> > >;

  /* normalizationParameters */
  using normParamROIMapType = std::map<std::string, scalarType>;   // ( 'min', v1),('max',v2),..
  using normParamType = std::map<std::string, normParamROIMapType>;

  /** set functions */
  void SetGradientSize( unsigned int length);

  void SetImagesOfInterestInOrder( WorkingImageVectorType& images);

  void SetImagesOfSpatialLocation(  std::map<std::string, WorkingImagePointer>& SpatialLocationImages);

  void SetCandidateROIs( std::map<std::string, WorkingImagePointer>& candidateROIMap);

  void SetROIInOrder( DataSet::StringVectorType roiInOrder );

  void SetFeatureInputOfROI( std::map<std::string, WorkingImagePointer>& featureImages );

  void SetInputVectorSize();

  unsigned int GetInputVectorSize();

  void SetLinearNormalizationOn(); // deprecated soon

  void SetNormalizationMethod( const std::string & normalizationMethod );

  void NormalizationOfVector( InputVectorMapType& currentFeatureVector, std::string ROIName );

  /** get function(s) */
  InputVectorMapType ComputeAndGetFeatureInputOfROI( std::string ROIName );

  /* HashGenerator From Index */
  /* hash function is based on fixed size of 'size'
   * This would not work if the size of image bigger than the size we are using here.
   * The size could be easily changed though.
   */
  static hashKeyType HashKeyFromIndex(const WorkingImageType::IndexType index);

  static WorkingImageType::IndexType HashIndexFromKey(const hashKeyType offSet);

private:
  int          m_gradientSize;
  unsigned int m_inputVectorSize;
  std::string  m_normalizationMethod;

  ImageLinearInterpolatorType::Pointer m_imageInterpolator;

  WorkingImageVectorType    m_imagesOfInterestInOrder;
  DataSet::StringVectorType m_roiIDsInOrder;

  /** deformed rho/phi/theta images*/
  std::map<std::string, WorkingImagePointer> m_spatialLocations;

  /** deformed candiateROIs */
  std::map<std::string, WorkingImagePointer> m_candidateROIs;

  /** gradient image of ROI */
  std::map<std::string, GradientImageType> m_gradientOfROI;

  /** feature output*/
  // std::map<std::string, InputVectorMapType> featureInputOfROI;

  /** m_normalization parameters*/
  /*  mapping from ROIname to the vector of mean/max
   *  for given serios of m_imagesOfInterestInOrder
   */
  std::map<std::string, minmaxPairVectorType> m_minmax;
  NormalizationParameterType                  m_statistics;

  std::map<std::string, FeatureNormalizationMethodEnum> m_featureNormalizationMap;

  /** private functions */
  // void ComputeFeatureInputOfROI( std::string ROIName);

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

  inline void AddFeaturesAlongGradient( std::string ROIName, const WorkingImagePointer& featureImage,
                                        WorkingImageType::IndexType currentPixelIndex,
                                        std::vector<scalarType>::iterator & elementIterator );

  inline std::map<std::string, scalarType> CalculateUnitDeltaAlongTheGradient(std::string ROIName,
                                                                              WorkingImageType::IndexType currentPixelIndex );

  inline std::pair<scalarType, scalarType>  SetMinMaxOfSubject( BinaryImageType::Pointer & labelImage,
                                                                const WorkingImagePointer & Image );

  inline float Sigmoid( const float x, const float t, const float r );

  inline float LinearScaling( float x, float min, float max );

  inline float doubleSigmoid( const float x, const float t, const float r1, const float r2);

  inline float ZScore( const float x, const float mu, const float sigma );
};
#endif
