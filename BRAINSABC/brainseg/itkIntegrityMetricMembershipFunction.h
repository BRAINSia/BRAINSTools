/*=========================================================================
 *
 *  Copyright Insight Software Consortium
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
/*
 * Author: Ali Ghayoor
 * at SINAPSE Lab,
 * The University of Iowa 2015
 */

#ifndef itkIntegrityMetricMembershipFunction_h
#define itkIntegrityMetricMembershipFunction_h

#include "itkVariableSizeMatrix.h"
#include "itkMeasurementVectorTraits.h"
#include "itkObject.h"
#include "itkObjectFactory.h"

namespace itk
{
namespace Statistics
{
/** \class IntegrityMetricMembershipFunction
 * \brief IntegrityMetricMembershipFunction models class membership
 * based on an Mahalanobis-weighted Euclidean distance.
 *
 * This class first calculates the Euclidean and the Mahalanobis distance
 * for each sample measurement from the computed mean point.
 * Then, Mahalanobis distances are scaled based on the maximum distance.
 * Scaled Mahalanobis distances are then used to wieght the calculated
 * Euclidean distances to create the output weighted distance vector.
 *
 * MahalanobisDistanceWeights = MahalanobisDistance/max(MahalanobisDistance)
 * OutputWeighteDistanceVector = EuclideanDistance * MahalanobisDistanceWeights
 *
 * The ``Evaluate()'' method returns true if all members of the OutputWeighteDistanceVector
 * are smaller than an input threshold value.
 *
 * \ingroup ITKStatistics
 */

template< typename TSample >
class IntegrityMetricMembershipFunction:
  public Object
{
public:
  /** Standard class typedefs */
  typedef IntegrityMetricMembershipFunction     Self;
  typedef Object                                Superclass;
  typedef SmartPointer< Self >                  Pointer;
  typedef SmartPointer< const Self >            ConstPointer;

  /** Strandard macros */
  itkTypeMacro(IntegrityMetricMembershipFunction, Object);
  itkNewMacro(Self);

  typedef TSample                                                     SampleType;

  /** Type of each measurement vector in sample */
  typedef typename SampleType::MeasurementVectorType                  MeasurementVectorType;

  /** Type of the length of each measurement vector */
  typedef typename SampleType::MeasurementVectorSizeType              MeasurementVectorSizeType;

  /** Type of measurement vector component value */
  typedef typename SampleType::MeasurementType                        MeasurementType;

  /** Type of a measurement vector, holding floating point values */
  typedef typename NumericTraits< MeasurementVectorType >::RealType   MeasurementVectorRealType;

  /** Type of a floating point measurement component value */
  typedef typename NumericTraits< MeasurementType >::RealType         MeasurementRealType;

  /** Type of the mean vector.  */
  typedef MeasurementVectorRealType                    MeanVectorType;

  /** Type of the covariance matrix */
  typedef VariableSizeMatrix< MeasurementRealType >    CovarianceMatrixType;

  /** Type of the output distance vector */
  typedef vnl_vector<double>                           DistanceVectorType;

  /** Set threshold */
  itkSetMacro(Threshold, float);

  /** Get the mean of the measurement samples. Mean is a vector type
   * similar to the measurement type but with a real element type. */
  itkGetConstReferenceMacro(Mean, MeanVectorType);

  /** Get the covariance matrix of the measurement samples. Covariance
   * matrix is a VariableSizeMatrix of real element type. */
  itkGetConstReferenceMacro(Covariance, CovarianceMatrixType);

  /**
   * Evaluate the weighted distance of a measurement using the
   * calculated Euclidean and Mahalanobis distances.
   * Note mean and covariance are computed to derive Mahalanobis distance.
   * This method returns true if all computed weighted distances are less
   * than input threshold. */
  bool Evaluate(const SampleType * measurementSample);

  /** Get calculated weighted distance vector */
  itkGetConstMacro(WeightedDistanceVector, DistanceVectorType);

  /** Set the length of the measurement vector. If this membership
   * function is templated over a vector type that can be resized,
   * the new size is set. If the vector type has a fixed size and an
   * attempt is made to change its size, an exception is
   * thrown. */
  virtual void SetMeasurementVectorSize(MeasurementVectorSizeType s)
  {
    // Test whether the vector type is resizable or not
  MeasurementVectorType m;

  if ( MeasurementVectorTraits::IsResizable(m) )
    {
    // then this is a resizable vector type
    //
    // if the new size is the same as the previou size, just return
    if ( s == this->m_MeasurementVectorSize )
      {
      return;
      }
    else
      {
      this->m_MeasurementVectorSize = s;
      this->Modified();
      }
    }
  else
    {
    // If this is a non-resizable vector type
    MeasurementVectorType     m3;
    MeasurementVectorSizeType defaultLength =
    NumericTraits<MeasurementVectorType>::GetLength(m3);
    // and the new length is different from the default one, then throw an exception
    if ( defaultLength != s )
      {
      itkExceptionMacro(
        "Attempting to change the measurement vector size of a non-resizable vector type" );
      }
    }
  }

  /** Get the length of the measurement vector */
  itkGetConstMacro(MeasurementVectorSize, MeasurementVectorSizeType);

protected:
  IntegrityMetricMembershipFunction();
  virtual ~IntegrityMetricMembershipFunction(void) {}
  void PrintSelf(std::ostream & os, Indent indent) const ITK_OVERRIDE;

  /** Set the mean used in the Mahalanobis distance.
    * This method run sanity checks after mean is computed.  */
  void SetMean(const MeanVectorType & mean);

  /** Set the covariance matrix.
    * This method run sanity checks after covariance is computed. */
  void SetCovariance(const CovarianceMatrixType & cov);

private:
  MeasurementVectorSizeType   m_MeasurementVectorSize;
  float                       m_Threshold;                // threshold value
  MeanVectorType              m_Mean;                     // mean
  CovarianceMatrixType        m_Covariance;               // covariance matrix
  DistanceVectorType          m_WeightedDistanceVector;   // output weighted distance vector
};
} // end of namespace Statistics
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkIntegrityMetricMembershipFunction.hxx"
#endif

#endif
