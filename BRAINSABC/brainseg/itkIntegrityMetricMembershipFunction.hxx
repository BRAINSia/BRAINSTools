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
#ifndef itkIntegrityMetricMembershipFunction_hxx
#define itkIntegrityMetricMembershipFunction_hxx

#include "itkIntegrityMetricMembershipFunction.h"

#include "itkVector.h"
#include "itkListSample.h"
#include "itkCovarianceSampleFilter.h"
#include "itkMahalanobisDistanceMetric.h"
#include "itkEuclideanDistanceMetric.h"

namespace itk
{
namespace Statistics
{
template <typename TSample>
IntegrityMetricMembershipFunction<TSample>::IntegrityMetricMembershipFunction()
{
  m_MeasurementVectorSize = NumericTraits<MeasurementVectorType>::GetLength(MeasurementVectorType());

  NumericTraits<MeanVectorType>::SetLength(m_Mean, this->GetMeasurementVectorSize());
  m_Mean.Fill(0.0f);

  m_Covariance.SetSize(this->GetMeasurementVectorSize(), this->GetMeasurementVectorSize());
  m_Covariance.SetIdentity();

  m_Threshold = 0.0;

  m_WeightedDistanceVector.set_size(0);
}

template <typename TSample>
void
IntegrityMetricMembershipFunction<TSample>::SetMean(const MeanVectorType & mean)
{
  if (this->GetMeasurementVectorSize())
  {
    MeasurementVectorTraits::Assert(mean,
                                    this->GetMeasurementVectorSize(),
                                    "IntegrityMetricMembershipFunction::SetMean(): Size of mean vector specified does "
                                    "not match the size of a measurement vector.");
  }
  else
  {
    // not already set, cache the size
    this->SetMeasurementVectorSize(mean.Size());
  }

  if (this->m_Mean != mean)
  {
    this->m_Mean = mean;
    this->Modified();
  }
}

template <typename TSample>
void
IntegrityMetricMembershipFunction<TSample>::SetCovariance(const CovarianceMatrixType & cov)
{
  // Sanity check
  if (cov.GetVnlMatrix().rows() != cov.GetVnlMatrix().cols())
  {
    itkExceptionMacro(<< "Covariance matrix must be square.");
  }

  if (this->GetMeasurementVectorSize())
  {
    if (cov.GetVnlMatrix().rows() != this->GetMeasurementVectorSize())
    {
      itkExceptionMacro(<< "Length of measurement vectors (" << this->GetMeasurementVectorSize()
                        << ") must be the same as the size of the covariance matrix (" << cov.GetVnlMatrix().rows()
                        << ").");
    }
  }
  else
  {
    // not already set, cache the size
    this->SetMeasurementVectorSize(cov.GetVnlMatrix().rows());
  }

  if (this->m_Covariance != cov)
  {
    this->m_Covariance = cov;
    this->Modified();
  }
}

template <typename TSample>
bool
IntegrityMetricMembershipFunction<TSample>::Evaluate(const SampleType * measurementSample)
{
  if (measurementSample->Size() == 0)
  {
    itkExceptionMacro(<< "There is no entry in input measurement sample.");
  }

  // If there is only one measurement sample per each plug area, that plug is pure!
  if (measurementSample->Size() == 1)
  {
    return true;
  }

  // compute mean and covariance
  using CovarianceAlgorithmType = itk::Statistics::CovarianceSampleFilter<SampleType>;
  typename CovarianceAlgorithmType::Pointer covarianceAlgorithm = CovarianceAlgorithmType::New();
  covarianceAlgorithm->SetInput(measurementSample);
  covarianceAlgorithm->Update();

  this->SetMean(covarianceAlgorithm->GetMean());
  this->SetCovariance(covarianceAlgorithm->GetCovarianceMatrix());

  // Compute Mahalanobis and Euclidean distances for each sample and put them in vectors
  using EDMetricType = itk::Statistics::EuclideanDistanceMetric<MeasurementVectorType>;
  typename EDMetricType::Pointer euclideanDist = EDMetricType::New();

  using MDMetricType = itk::Statistics::MahalanobisDistanceMetric<MeasurementVectorType>;
  typename MDMetricType::Pointer mahalanobisDist = MDMetricType::New();
  mahalanobisDist->SetMeasurementVectorSize(this->GetMeasurementVectorSize());
  mahalanobisDist->SetCovariance(this->GetCovariance().GetVnlMatrix());

  vnl_vector<double> ed_vector(measurementSample->Size()); // vector including euclidean distances for each sample
  vnl_vector<double> md_vector(measurementSample->Size()); // vector including mahalanobis distances for each sample
  this->m_WeightedDistanceVector.clear();

  unsigned int i = 0;
  for (typename SampleType::ConstIterator s_iter = measurementSample->Begin(); s_iter != measurementSample->End();
       ++s_iter)
  {
    ed_vector[i] = euclideanDist->Evaluate(this->GetMean(), s_iter.GetMeasurementVector());
    md_vector[i] = mahalanobisDist->Evaluate(this->GetMean(), s_iter.GetMeasurementVector());
    ++i;
  }

  // If points are so close to the mean, you can rest assure they belong to a pure plug.
  // However, this ideal case mostly happen in background area!
  if (ed_vector.max_value() < 1e-2)
  {
    return true;
  }

  md_vector /= md_vector.max_value(); // Normalize mahalanobis distances to the maximum distance
  this->m_WeightedDistanceVector = element_product(ed_vector, md_vector);

  bool belongs = (this->m_WeightedDistanceVector.max_value() < this->m_Threshold);
  return belongs;
}

template <typename TSample>
void
IntegrityMetricMembershipFunction<TSample>::PrintSelf(std::ostream & os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);
  os << indent << "Measurement vector size: " << this->m_MeasurementVectorSize << std::endl;
  os << indent << "Mean: " << this->m_Mean << std::endl;
  os << indent << "Covariance: " << std::endl;
  os << this->m_Covariance.GetVnlMatrix();
  os << indent << "Weighted distance matrix: " << std::endl;
  os << this->m_WeightedDistanceVector << std::endl;
}

} // end namespace Statistics
} // end of namespace itk

#endif
