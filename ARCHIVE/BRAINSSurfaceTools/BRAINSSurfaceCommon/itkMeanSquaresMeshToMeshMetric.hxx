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
#ifndef __itkMeanSquaresMeshToMeshMetric_hxx
#define __itkMeanSquaresMeshToMeshMetric_hxx

// First make sure that the configuration is available.
// This line can be removed once the optimized versions
// gets integrated into the main directories.
#include "itkConfigure.h"

#include "itkMeanSquaresMeshToMeshMetric.h"

namespace itk
{
/**
 * Constructor
 */
template <typename TFixedMesh, typename TMovingMesh>
MeanSquaresMeshToMeshMetric<TFixedMesh, TMovingMesh>::MeanSquaresMeshToMeshMetric()
{
  itkDebugMacro("Constructor");
}

/**
 * Get the match Measure
 */
template <typename TFixedMesh, typename TMovingMesh>
typename MeanSquaresMeshToMeshMetric<TFixedMesh, TMovingMesh>::MeasureType
MeanSquaresMeshToMeshMetric<TFixedMesh, TMovingMesh>::GetValue(const TransformParametersType & parameters) const
{
  FixedMeshConstPointer fixedMesh = this->GetFixedMesh();

  if (!fixedMesh)
  {
    itkExceptionMacro(<< "Fixed point set has not been assigned");
  }

  PointIterator pointItr = fixedMesh->GetPoints()->Begin();
  PointIterator pointEnd = fixedMesh->GetPoints()->End();

  PointDataIterator pointDataItr = fixedMesh->GetPointData()->Begin();
  PointDataIterator pointDataEnd = fixedMesh->GetPointData()->End();

  MeasureType sumOfSquaresDifferences = NumericTraits<MeasureType>::ZeroValue();

  this->m_NumberOfPixelsCounted = 0;

  this->SetTransformParameters(parameters);

  using InterpolationPointType = typename InterpolatorType::PointType;
  InterpolationPointType pointToEvaluate;

  while (pointItr != pointEnd && pointDataItr != pointDataEnd)
  {
    InputPointType inputPoint;
    inputPoint.CastFrom(pointItr.Value());

    if (this->m_FixedMask && !this->m_FixedMask->IsInsideInWorldSpace(inputPoint))
    {
      ++pointItr;
      ++pointDataItr;
      continue;
    }

    OutputPointType transformedPoint = this->m_Transform->TransformPoint(inputPoint);

    if (this->m_MovingMask && !this->m_MovingMask->IsInsideInWorldSpace(transformedPoint))
    {
      ++pointItr;
      ++pointDataItr;
      continue;
    }

    // FIXME: if( this->m_Interpolator->IsInsideSurface( transformedPoint ) )
    {
      pointToEvaluate.CastFrom(transformedPoint);
      const RealDataType movingValue = this->m_Interpolator->Evaluate(pointToEvaluate);
      const RealDataType fixedValue = pointDataItr.Value();
      const RealDataType diff = movingValue - fixedValue;
      sumOfSquaresDifferences += diff * diff;
      this->m_NumberOfPixelsCounted++;
    }

    ++pointItr;
    ++pointDataItr;
  }

  if (!this->m_NumberOfPixelsCounted)
  {
    itkExceptionMacro(<< "All the points mapped to outside of the moving image");
  }

  const double averageOfSquaredDifferences = sumOfSquaresDifferences / this->m_NumberOfPixelsCounted;

  return averageOfSquaredDifferences;
}

/**
 * Get the Derivative Measure
 */
template <typename TFixedMesh, typename TMovingMesh>
void
MeanSquaresMeshToMeshMetric<TFixedMesh, TMovingMesh>::GetDerivative(const TransformParametersType & parameters,
                                                                    DerivativeType &                derivative) const
{
  itkDebugMacro("GetDerivative( " << parameters << " ) ");

  FixedMeshConstPointer fixedMesh = this->m_FixedMesh;

  if (!fixedMesh)
  {
    itkExceptionMacro(<< "Fixed image has not been assigned");
  }

  PointIterator pointItr = fixedMesh->GetPoints()->Begin();
  PointIterator pointEnd = fixedMesh->GetPoints()->End();

  PointDataIterator pointDataItr = fixedMesh->GetPointData()->Begin();
  PointDataIterator pointDataEnd = fixedMesh->GetPointData()->End();

  this->m_NumberOfPixelsCounted = 0;

  this->SetTransformParameters(parameters);

  const unsigned int ParametersDimension = this->GetNumberOfParameters();
  derivative = DerivativeType(ParametersDimension);
  derivative.Fill(NumericTraits<typename DerivativeType::ValueType>::ZeroValue());

  using InterpolationPointType = typename InterpolatorType::PointType;
  InterpolationPointType pointToEvaluate;

  while (pointItr != pointEnd && pointDataItr != pointDataEnd)
  {
    InputPointType inputPoint;
    inputPoint.CastFrom(pointItr.Value());

    if (this->m_FixedMask && !this->m_FixedMask->IsInsideInWorldSpace(inputPoint))
    {
      ++pointItr;
      ++pointDataItr;
      continue;
    }

    OutputPointType transformedPoint = this->m_Transform->TransformPoint(inputPoint);

    if (this->m_MovingMask && !this->m_MovingMask->IsInsideInWorldSpace(transformedPoint))
    {
      ++pointItr;
      ++pointDataItr;
      continue;
    }

    // FIXME:  if( this->m_Interpolator->IsInsideBuffer( transformedPoint ) )
    {
      pointToEvaluate.CastFrom(transformedPoint);
      const RealDataType movingValue = this->m_Interpolator->Evaluate(pointToEvaluate);

      TransformJacobianType jacobian;
      this->m_Transform->ComputeJacobianWithRespectToParameters(inputPoint, jacobian);

      const RealDataType fixedValue = pointDataItr.Value();
      this->m_NumberOfPixelsCounted++;

      const RealDataType diff = movingValue - fixedValue;

      DerivativeDataType gradient;
      this->m_Interpolator->EvaluateDerivative(pointToEvaluate, gradient);
      for (unsigned int par = 0; par < ParametersDimension; par++)
      {
        RealDataType sum = NumericTraits<RealDataType>::ZeroValue();
        for (unsigned int dim = 0; dim < MovingMeshDimension; dim++)
        {
          sum += 2.0 * diff * jacobian(dim, par) * gradient[dim];
        }
        derivative[par] += sum;
      }
    }

    ++pointItr;
    ++pointDataItr;
  }

  if (!this->m_NumberOfPixelsCounted)
  {
    itkExceptionMacro(<< "All the points mapped to outside of the moving image");
  }
  else
  {
    for (unsigned int i = 0; i < ParametersDimension; i++)
    {
      derivative[i] /= this->m_NumberOfPixelsCounted;
    }
  }
}

/*
 * Get both the match Measure and theDerivative Measure
 */
template <typename TFixedMesh, typename TMovingMesh>
void
MeanSquaresMeshToMeshMetric<TFixedMesh, TMovingMesh>::GetValueAndDerivative(const TransformParametersType & parameters,
                                                                            MeasureType &                   value,
                                                                            DerivativeType & derivative) const
{
  itkDebugMacro("GetValueAndDerivative( " << parameters << " ) ");

  FixedMeshConstPointer fixedMesh = this->m_FixedMesh;

  if (!fixedMesh)
  {
    itkExceptionMacro(<< "Fixed image has not been assigned");
  }

  PointIterator pointItr = fixedMesh->GetPoints()->Begin();
  PointIterator pointEnd = fixedMesh->GetPoints()->End();

  PointDataIterator pointDataItr = fixedMesh->GetPointData()->Begin();
  PointDataIterator pointDataEnd = fixedMesh->GetPointData()->End();

  MeasureType sumOfSquaresDifferences = NumericTraits<MeasureType>::ZeroValue();

  this->m_NumberOfPixelsCounted = 0;

  this->SetTransformParameters(parameters);

  const unsigned int ParametersDimension = this->GetNumberOfParameters();
  derivative = DerivativeType(ParametersDimension);
  derivative.Fill(NumericTraits<typename DerivativeType::ValueType>::ZeroValue());

  using InterpolationPointType = typename InterpolatorType::PointType;
  InterpolationPointType pointToEvaluate;

  while (pointItr != pointEnd && pointDataItr != pointDataEnd)
  {
    InputPointType inputPoint;
    inputPoint.CastFrom(pointItr.Value());

    if (this->m_FixedMask && !this->m_FixedMask->IsInsideInWorldSpace(inputPoint))
    {
      ++pointItr;
      ++pointDataItr;
      continue;
    }

    OutputPointType transformedPoint = this->m_Transform->TransformPoint(inputPoint);

    if (this->m_MovingMask && !this->m_MovingMask->IsInsideInWorldSpace(transformedPoint))
    {
      ++pointItr;
      ++pointDataItr;
      continue;
    }

    // FIXME:  if( this->m_Interpolator->IsInsideBuffer( transformedPoint ) )
    {
      pointToEvaluate.CastFrom(transformedPoint);
      const RealDataType movingValue = this->m_Interpolator->Evaluate(pointToEvaluate);

      TransformJacobianType jacobian;
      this->m_Transform->ComputeJacobianWithRespectToParameters(inputPoint, jacobian);

      const RealDataType fixedValue = pointDataItr.Value();
      this->m_NumberOfPixelsCounted++;

      const RealDataType diff = movingValue - fixedValue;

      sumOfSquaresDifferences += diff * diff;

      DerivativeDataType gradient;

      this->m_Interpolator->EvaluateDerivative(pointToEvaluate, gradient);
      for (unsigned int par = 0; par < ParametersDimension; par++)
      {
        RealDataType sum = NumericTraits<RealDataType>::ZeroValue();
        for (unsigned int dim = 0; dim < MovingMeshDimension; dim++)
        {
          sum += 2.0 * diff * jacobian(dim, par) * gradient[dim];
        }
        derivative[par] += sum;
      }
    }

    ++pointItr;
    ++pointDataItr;
  }

  if (!this->m_NumberOfPixelsCounted)
  {
    itkExceptionMacro(<< "All the points mapped to outside of the moving image");
  }
  for (unsigned int i = 0; i < ParametersDimension; i++)
  {
    derivative[i] /= this->m_NumberOfPixelsCounted;
  }
  const double averageOfSquaredDifferences = sumOfSquaresDifferences / this->m_NumberOfPixelsCounted;

  value = averageOfSquaredDifferences;
}
} // end namespace itk

#endif
