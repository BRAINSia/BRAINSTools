/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkMeanSquaresMeshToMeshMetric.txx,v $
  Language:  C++
  Date:      $Date:  $
  Version:   $Revision:  $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
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
template <class TFixedMesh, class TMovingMesh>
MeanSquaresMeshToMeshMetric<TFixedMesh, TMovingMesh>
::MeanSquaresMeshToMeshMetric()
{
  itkDebugMacro("Constructor");
}

/**
 * Get the match Measure
 */
template <class TFixedMesh, class TMovingMesh>
typename MeanSquaresMeshToMeshMetric<TFixedMesh, TMovingMesh>::MeasureType
MeanSquaresMeshToMeshMetric<TFixedMesh, TMovingMesh>
::GetValue( const TransformParametersType & parameters ) const
{
  FixedMeshConstPointer fixedMesh = this->GetFixedMesh();

  if( !fixedMesh )
    {
    itkExceptionMacro( << "Fixed point set has not been assigned" );
    }

  PointIterator pointItr = fixedMesh->GetPoints()->Begin();
  PointIterator pointEnd = fixedMesh->GetPoints()->End();

  PointDataIterator pointDataItr = fixedMesh->GetPointData()->Begin();
  PointDataIterator pointDataEnd = fixedMesh->GetPointData()->End();

  MeasureType sumOfSquaresDifferences = NumericTraits<MeasureType>::Zero;

  this->m_NumberOfPixelsCounted = 0;

  this->SetTransformParameters( parameters );

  typedef typename InterpolatorType::PointType InterpolationPointType;
  InterpolationPointType pointToEvaluate;

  while( pointItr != pointEnd && pointDataItr != pointDataEnd )
    {
    InputPointType inputPoint;
    inputPoint.CastFrom( pointItr.Value() );

    if( this->m_FixedMask && !this->m_FixedMask->IsInside( inputPoint ) )
      {
      ++pointItr;
      ++pointDataItr;
      continue;
      }

    OutputPointType transformedPoint = this->m_Transform->TransformPoint( inputPoint );

    if( this->m_MovingMask && !this->m_MovingMask->IsInside( transformedPoint ) )
      {
      ++pointItr;
      ++pointDataItr;
      continue;
      }

    // FIXME: if( this->m_Interpolator->IsInsideSurface( transformedPoint ) )
      {
      pointToEvaluate.CastFrom( transformedPoint );
      const RealDataType movingValue  = this->m_Interpolator->Evaluate( pointToEvaluate );
      const RealDataType fixedValue   = pointDataItr.Value();
      const RealDataType diff = movingValue - fixedValue;
      sumOfSquaresDifferences += diff * diff;
      this->m_NumberOfPixelsCounted++;
      }

    ++pointItr;
    ++pointDataItr;
    }

  if( !this->m_NumberOfPixelsCounted )
    {
    itkExceptionMacro(<< "All the points mapped to outside of the moving image");
    }

  const double averageOfSquaredDifferences = sumOfSquaresDifferences / this->m_NumberOfPixelsCounted;

  return averageOfSquaredDifferences;
}

/**
 * Get the Derivative Measure
 */
template <class TFixedMesh, class TMovingMesh>
void
MeanSquaresMeshToMeshMetric<TFixedMesh, TMovingMesh>
::GetDerivative( const TransformParametersType & parameters,
                 DerivativeType & derivative  ) const
{
  itkDebugMacro("GetDerivative( " << parameters << " ) ");

  FixedMeshConstPointer fixedMesh = this->m_FixedMesh;

  if( !fixedMesh )
    {
    itkExceptionMacro( << "Fixed image has not been assigned" );
    }

  PointIterator pointItr = fixedMesh->GetPoints()->Begin();
  PointIterator pointEnd = fixedMesh->GetPoints()->End();

  PointDataIterator pointDataItr = fixedMesh->GetPointData()->Begin();
  PointDataIterator pointDataEnd = fixedMesh->GetPointData()->End();

  this->m_NumberOfPixelsCounted = 0;

  this->SetTransformParameters( parameters );

  const unsigned int ParametersDimension = this->GetNumberOfParameters();
  derivative = DerivativeType( ParametersDimension );
  derivative.Fill( NumericTraits<ITK_TYPENAME DerivativeType::ValueType>::Zero );

  typedef typename InterpolatorType::PointType InterpolationPointType;
  InterpolationPointType pointToEvaluate;

  while( pointItr != pointEnd && pointDataItr != pointDataEnd )
    {
    InputPointType inputPoint;
    inputPoint.CastFrom( pointItr.Value() );

    if( this->m_FixedMask && !this->m_FixedMask->IsInside( inputPoint ) )
      {
      ++pointItr;
      ++pointDataItr;
      continue;
      }

    OutputPointType transformedPoint = this->m_Transform->TransformPoint( inputPoint );

    if( this->m_MovingMask && !this->m_MovingMask->IsInside( transformedPoint ) )
      {
      ++pointItr;
      ++pointDataItr;
      continue;
      }

    // FIXME:  if( this->m_Interpolator->IsInsideBuffer( transformedPoint ) )
      {
      pointToEvaluate.CastFrom( transformedPoint );
      const RealDataType movingValue  = this->m_Interpolator->Evaluate( pointToEvaluate );

      const TransformJacobianType & jacobian =
        this->m_Transform->GetJacobian( inputPoint );

      const RealDataType fixedValue     = pointDataItr.Value();
      this->m_NumberOfPixelsCounted++;

      const RealDataType diff = movingValue - fixedValue;

      DerivativeDataType gradient;

      this->m_Interpolator->EvaluateDerivative( pointToEvaluate, gradient );
      for( unsigned int par = 0; par < ParametersDimension; par++ )
        {
        RealDataType sum = NumericTraits<RealDataType>::Zero;
        for( unsigned int dim = 0; dim < MovingMeshDimension; dim++ )
          {
          sum += 2.0 * diff * jacobian( dim, par ) * gradient[dim];
          }
        derivative[par] += sum;
        }
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
    }
}

/*
 * Get both the match Measure and theDerivative Measure
 */
template <class TFixedMesh, class TMovingMesh>
void
MeanSquaresMeshToMeshMetric<TFixedMesh, TMovingMesh>
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

  MeasureType sumOfSquaresDifferences = NumericTraits<MeasureType>::Zero;

  this->m_NumberOfPixelsCounted = 0;

  this->SetTransformParameters( parameters );

  const unsigned int ParametersDimension = this->GetNumberOfParameters();
  derivative = DerivativeType( ParametersDimension );
  derivative.Fill( NumericTraits<ITK_TYPENAME DerivativeType::ValueType>::Zero );

  typedef typename InterpolatorType::PointType InterpolationPointType;
  InterpolationPointType pointToEvaluate;

  while( pointItr != pointEnd && pointDataItr != pointDataEnd )
    {
    InputPointType inputPoint;
    inputPoint.CastFrom( pointItr.Value() );

    if( this->m_FixedMask && !this->m_FixedMask->IsInside( inputPoint ) )
      {
      ++pointItr;
      ++pointDataItr;
      continue;
      }

    OutputPointType transformedPoint = this->m_Transform->TransformPoint( inputPoint );

    if( this->m_MovingMask && !this->m_MovingMask->IsInside( transformedPoint ) )
      {
      ++pointItr;
      ++pointDataItr;
      continue;
      }

    // FIXME:  if( this->m_Interpolator->IsInsideBuffer( transformedPoint ) )
      {
      pointToEvaluate.CastFrom( transformedPoint );
      const RealDataType movingValue  = this->m_Interpolator->Evaluate( pointToEvaluate );

      const TransformJacobianType & jacobian =
        this->m_Transform->GetJacobian( inputPoint );

      const RealDataType fixedValue     = pointDataItr.Value();
      this->m_NumberOfPixelsCounted++;

      const RealDataType diff = movingValue - fixedValue;

      sumOfSquaresDifferences += diff * diff;

      DerivativeDataType gradient;

      this->m_Interpolator->EvaluateDerivative( pointToEvaluate, gradient );
      for( unsigned int par = 0; par < ParametersDimension; par++ )
        {
        RealDataType sum = NumericTraits<RealDataType>::Zero;
        for( unsigned int dim = 0; dim < MovingMeshDimension; dim++ )
          {
          sum += 2.0 * diff * jacobian( dim, par ) * gradient[dim];
          }
        derivative[par] += sum;
        }
      }

    ++pointItr;
    ++pointDataItr;
    }

  if( !this->m_NumberOfPixelsCounted )
    {
    itkExceptionMacro(<< "All the points mapped to outside of the moving image");
    }
  for( unsigned int i = 0; i < ParametersDimension; i++ )
    {
    derivative[i] /= this->m_NumberOfPixelsCounted;
    }
  const double averageOfSquaredDifferences = sumOfSquaresDifferences / this->m_NumberOfPixelsCounted;

  value = averageOfSquaredDifferences;
}
} // end namespace itk

#endif
