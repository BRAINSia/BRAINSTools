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
#ifndef __itkAnalyticalMeshToMeshMetric_hxx
#define __itkAnalyticalMeshToMeshMetric_hxx

// First make sure that the configuration is available.
// This line can be removed once the optimized versions
// gets integrated into the main directories.
#include "itkConfigure.h"

#include "itkAnalyticalMeshToMeshMetric.h"

namespace itk
{
/**
 * Constructor
 */
template <class TFixedMesh, class TMovingMesh>
AnalyticalMeshToMeshMetric<TFixedMesh, TMovingMesh>
::AnalyticalMeshToMeshMetric()
{
  itkDebugMacro("Constructor");
}

/**
 * Get the match Measure
 */
template <class TFixedMesh, class TMovingMesh>
typename AnalyticalMeshToMeshMetric<TFixedMesh, TMovingMesh>::MeasureType
AnalyticalMeshToMeshMetric<TFixedMesh, TMovingMesh>
::GetValue( const TransformParametersType & parameters ) const
{
  // here the value of the "metric" is simply the moving value...
  FixedMeshConstPointer fixedMesh = this->GetFixedMesh();

  if( !fixedMesh )
    {
    itkExceptionMacro( << "Fixed point set has not been assigned" );
    }

  PointIterator pointItr = fixedMesh->GetPoints()->Begin();
  PointIterator pointEnd = fixedMesh->GetPoints()->End();

  PointDataIterator pointDataItr = fixedMesh->GetPointData()->Begin();
  PointDataIterator pointDataEnd = fixedMesh->GetPointData()->End();

  MeasureType measure = NumericTraits<MeasureType>::Zero;

  this->m_NumberOfPixelsCounted = 0;

  this->SetTransformParameters( parameters );

  while( pointItr != pointEnd && pointDataItr != pointDataEnd )
    {
    InputPointType inputPoint;
    inputPoint.CastFrom( pointItr.Value() );
    OutputPointType transformedPoint =
      this->m_Transform->TransformPoint( inputPoint );

    // FIXME: if( this->m_Interpolator->IsInsideSurface( transformedPoint ) )
      {
      const RealDataType movingValue  = this->m_Interpolator->Evaluate( transformedPoint );
      const RealDataType fixedValue   = pointDataItr.Value();
      const RealDataType diff = movingValue - fixedValue;
      measure += diff * diff;
      this->m_NumberOfPixelsCounted++;
      }

    ++pointItr;
    ++pointDataItr;
    }

  if( !this->m_NumberOfPixelsCounted )
    {
    itkExceptionMacro(<< "All the points mapped to outside of the moving image");
    }
  else
    {
    measure /= this->m_NumberOfPixelsCounted;
    }

  return measure;
}

/**
 * Get the Derivative Measure
 */
template <class TFixedMesh, class TMovingMesh>
void
AnalyticalMeshToMeshMetric<TFixedMesh, TMovingMesh>
::GetDerivative( const TransformParametersType & parameters,
                 DerivativeType & derivative  ) const
{
#ifdef MICHEL_LATER

  itkDebugMacro("GetDerivative( " << parameters << " ) ");

  if( !this->GetGradientMesh() )
    {
    itkExceptionMacro(<< "The gradient image is null, maybe you forgot to call Initialize()");
    }

  FixedMeshConstPointer fixedMesh = this->m_FixedMesh;

  if( !fixedMesh )
    {
    itkExceptionMacro( << "Fixed image has not been assigned" );
    }

  const unsigned int MeshDimension = FixedMeshType::MeshDimension;

  typedef  itk::MeshRegionConstIteratorWithIndex<
      FixedMeshType> FixedIteratorType;

  FixedIteratorType ti( fixedMesh, this->GetFixedMeshRegion() );

  typename FixedMeshType::IndexType index;

  this->m_NumberOfPixelsCounted = 0;

  this->SetTransformParameters( parameters );

  const unsigned int ParametersDimension = this->GetNumberOfParameters();
  derivative = DerivativeType( ParametersDimension );
  derivative.Fill( NumericTraits<ITK_TYPENAME DerivativeType::ValueType>::Zero );

  ti.GoToBegin();

  while( !ti.IsAtEnd() )
    {
    index = ti.GetIndex();

    InputPointType inputPoint;
    fixedMesh->TransformIndexToPhysicalPoint( index, inputPoint );

    if( this->m_FixedMeshMask && !this->m_FixedMeshMask->IsInside( inputPoint ) )
      {
      ++ti;
      continue;
      }

    OutputPointType transformedPoint = this->m_Transform->TransformPoint( inputPoint );

    if( this->m_MovingMeshMask && !this->m_MovingMeshMask->IsInside( transformedPoint ) )
      {
      ++ti;
      continue;
      }

    if( this->m_Interpolator->IsInsideBuffer( transformedPoint ) )
      {
      const RealDataType movingValue  = this->m_Interpolator->Evaluate( transformedPoint );

      const TransformJacobianType & jacobian =
        this->m_Transform->GetJacobian( inputPoint );

      const RealDataType fixedValue     = ti.Value();
      this->m_NumberOfPixelsCounted++;
      const RealDataType diff = movingValue - fixedValue;

      // Get the gradient by NearestNeighboorInterpolation:
      // which is equivalent to round up the point components.
      typedef typename OutputPointType::CoordRepType CoordRepType;
      typedef ContinuousIndex<CoordRepType, MovingMeshType::MeshDimension>
        MovingMeshContinuousIndexType;

      MovingMeshContinuousIndexType tempIndex;
      this->m_MovingMesh->TransformPhysicalPointToContinuousIndex( transformedPoint, tempIndex );

      typename MovingMeshType::IndexType mappedIndex;
      mappedIndex.CopyWithRound( tempIndex );

      DerivativeType gradient;

      this->m_Interpolator->EvaluateDerivative( transformedPoint, gradient );
      for( unsigned int par = 0; par < ParametersDimension; par++ )
        {
        RealDataType sum = NumericTraits<RealDataType>::Zero;
        for( unsigned int dim = 0; dim < MeshDimension; dim++ )
          {
          sum += 2.0 * diff * jacobian( dim, par ) * gradient[dim];
          }
        derivative[par] += sum;
        }
      }

    ++ti;
    }

  if( !this->m_NumberOfPixelsCounted )
    {
    itkExceptionMacro(<< "All the points mapped to outside of the moving image");
    }
  else
    {
    for( unsigned int i = 0; i < ParametersDimension; i++ )
      {
      derivative[i] /= this->m_NumberOfPixelsCounted;
      }
    }

#endif
}

/*
 * Get both the match Measure and theDerivative Measure
 */
template <class TFixedMesh, class TMovingMesh>
void
AnalyticalMeshToMeshMetric<TFixedMesh, TMovingMesh>
::GetValueAndDerivative(const TransformParametersType & parameters,
                        MeasureType & value, DerivativeType  & derivative) const
{
  itkDebugMacro("GetValueAndDerivative( " << parameters << " ) ");

  FixedMeshConstPointer fixedMesh = this->m_FixedMesh;

  if( !fixedMesh )
    {
    itkExceptionMacro( << "Fixed image has not been assigned" );
    }

  PointIterator pointItr = fixedMesh->GetPoints()->Begin();
  PointIterator pointEnd = fixedMesh->GetPoints()->End();

  PointDataIterator pointDataItr = fixedMesh->GetPointData()->Begin();
  PointDataIterator pointDataEnd = fixedMesh->GetPointData()->End();

  MeasureType measure = NumericTraits<MeasureType>::Zero;

  this->m_NumberOfPixelsCounted = 0;

  this->SetTransformParameters( parameters );

  const unsigned int ParametersDimension = this->GetNumberOfParameters();
  derivative = DerivativeType( ParametersDimension );
  derivative.Fill( NumericTraits<ITK_TYPENAME DerivativeType::ValueType>::Zero );

  while( pointItr != pointEnd && pointDataItr != pointDataEnd )
    {
    InputPointType inputPoint;
    inputPoint.CastFrom( pointItr.Value() );

    if( this->m_FixedMeshMask && !this->m_FixedMeshMask->IsInside( inputPoint ) )
      {
      ++pointItr;
      ++pointDataItr;
      continue;
      }

    OutputPointType transformedPoint = this->m_Transform->TransformPoint( inputPoint );

    if( this->m_MovingMeshMask && !this->m_MovingMeshMask->IsInside( transformedPoint ) )
      {
      ++pointItr;
      ++pointDataItr;
      continue;
      }

    // FIXME:  if( this->m_Interpolator->IsInsideBuffer( transformedPoint ) )
      {
      const RealDataType movingValue  = this->m_Interpolator->Evaluate( transformedPoint );

      const TransformJacobianType & jacobian =
        this->m_Transform->GetJacobian( inputPoint );

      const RealDataType fixedValue     = pointDataItr.Value();
      this->m_NumberOfPixelsCounted++;

      const RealDataType diff = movingValue - fixedValue;

      measure += diff * diff;

      DerivativeDataType gradient;

      this->m_Interpolator->EvaluateDerivative( transformedPoint, gradient );
      for( unsigned int par = 0; par < ParametersDimension; par++ )
        {
        RealDataType sum = NumericTraits<RealDataType>::Zero;
        for( unsigned int dim = 0; dim < MovingMeshDimension; dim++ )
          {
          sum += 2.0 * diff * jacobian( dim, par ) * gradient[dim];
          }
        derivative[par] += sum;
        }

#ifdef MICHELLATER
      std::cout << "  Transform " << this->m_Transform->GetParameters() << std::endl;
      std::cout << "  inputPoint " << inputPoint << std::endl;
      std::cout << "  transformedPoint " << transformedPoint << std::endl;
      std::cout << "  gradient " << gradient  << std::endl;
      std::cout << "  jacobian " << jacobian << std::endl;
      std::cout << "  derivative " << derivative << std::endl;
      std::cout << std::endl;
#endif
      }

    ++pointItr;
    ++pointDataItr;
    }

  if( !this->m_NumberOfPixelsCounted )
    {
    itkExceptionMacro(<< "All the points mapped to outside of the moving image");
    }
  else
    {
    for( unsigned int i = 0; i < ParametersDimension; i++ )
      {
      derivative[i] /= this->m_NumberOfPixelsCounted;
      }
    measure /= this->m_NumberOfPixelsCounted;
    }

  value = measure;
}
} // end namespace itk

#endif
