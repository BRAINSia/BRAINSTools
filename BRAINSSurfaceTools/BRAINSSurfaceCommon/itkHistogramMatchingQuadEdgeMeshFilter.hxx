/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkHistogramMatchingQuadEdgeMeshFilter.h,v $
  Language:  C++
  Date:      $Date: 2010-05-14 09:26:05 $
  Version:   $Revision: 1.21 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkHistogramMatchingQuadEdgeMeshFilter_hxx
#define __itkHistogramMatchingQuadEdgeMeshFilter_hxx

#include "itkHistogramMatchingQuadEdgeMeshFilter.h"
#include "itkNumericTraits.h"
#include <vector>

namespace itk
{
template <class TInputMesh, class TOutputMesh, class THistogramMeasurement>
HistogramMatchingQuadEdgeMeshFilter<TInputMesh, TOutputMesh, THistogramMeasurement>
::HistogramMatchingQuadEdgeMeshFilter()
{
  this->SetNumberOfRequiredInputs( 2 );
  this->SetNumberOfRequiredOutputs( 1 );
  this->SetNumberOfIndexedOutputs( 1 );

  this->SetNthOutput( 0, OutputMeshType::New() );

  m_NumberOfHistogramLevels = 256;
  m_NumberOfMatchPoints = 1;

  m_QuantileTable.set_size( 3, m_NumberOfMatchPoints + 2 );
  m_QuantileTable.fill(0);
  m_Gradients.set_size( m_NumberOfMatchPoints + 1 );
  m_Gradients.fill(0);

  m_LowerGradient = 0.0;
  m_UpperGradient = 0.0;

  // Create histograms.
  m_SourceHistogram = HistogramType::New();
  m_ReferenceHistogram = HistogramType::New();
  m_OutputHistogram = HistogramType::New();
}

template <class TInputMesh, class TOutputMesh, class THistogramMeasurement>
HistogramMatchingQuadEdgeMeshFilter<TInputMesh, TOutputMesh, THistogramMeasurement>
::~HistogramMatchingQuadEdgeMeshFilter()
{
}

template <class TInputMesh, class TOutputMesh, class THistogramMeasurement>
void
HistogramMatchingQuadEdgeMeshFilter<TInputMesh, TOutputMesh, THistogramMeasurement>
::SetReferenceMesh( const InputMeshType * reference )
{
  itkDebugMacro("setting reference Mesh to " << reference);
  if( reference != static_cast<const InputMeshType *>(this->ProcessObject::GetInput( 1 ) ) )
    {
    this->ProcessObject::SetNthInput(1, const_cast<InputMeshType *>( reference ) );
    this->Modified();
    }
}

template <class TInputMesh, class TOutputMesh, class THistogramMeasurement>
const typename
HistogramMatchingQuadEdgeMeshFilter<TInputMesh, TOutputMesh, THistogramMeasurement>::InputMeshType
* HistogramMatchingQuadEdgeMeshFilter<TInputMesh, TOutputMesh, THistogramMeasurement>
::GetReferenceMesh() const
  {
  Self *                surrogate = const_cast<Self *>( this );
  const InputMeshType * reference =
    static_cast<const InputMeshType *>( surrogate->ProcessObject::GetInput(1) );
  return reference;
  }

template <class TInputMesh, class TOutputMesh, class THistogramMeasurement>
void
HistogramMatchingQuadEdgeMeshFilter<TInputMesh, TOutputMesh, THistogramMeasurement>
::SetSourceMesh( const InputMeshType * source )
{
  this->SetInput( source );
}

template <class TInputMesh, class TOutputMesh, class THistogramMeasurement>
const typename HistogramMatchingQuadEdgeMeshFilter<TInputMesh, TOutputMesh, THistogramMeasurement>::InputMeshType
* HistogramMatchingQuadEdgeMeshFilter<TInputMesh, TOutputMesh, THistogramMeasurement>
::GetSourceMesh( void ) const
  {
  return this->GetInput();
  }

template <class TInputMesh, class TOutputMesh, class THistogramMeasurement>
void
HistogramMatchingQuadEdgeMeshFilter<TInputMesh, TOutputMesh, THistogramMeasurement>
::BeforeTransform()
{
  this->CopyInputMeshToOutputMesh();

  InputMeshConstPointer source    = this->GetSourceMesh();
  InputMeshConstPointer reference = this->GetReferenceMesh();

  this->ComputeMinMax( source, m_SourceMinValue,
                       m_SourceMaxValue );
  this->ComputeMinMax( reference, m_ReferenceMinValue,
                       m_ReferenceMaxValue );

  this->ConstructHistogram( source, m_SourceHistogram,
                            m_SourceMinValue, m_SourceMaxValue );
  this->ConstructHistogram( reference, m_ReferenceHistogram,
                            m_ReferenceMinValue,
                            m_ReferenceMaxValue );

  // Fill in the quantile table.
  m_QuantileTable.set_size( 3, m_NumberOfMatchPoints + 2 );
  m_QuantileTable[0][0] = m_SourceMinValue;
  m_QuantileTable[1][0] = m_ReferenceMinValue;

  m_QuantileTable[0][m_NumberOfMatchPoints + 1] = m_SourceMaxValue;
  m_QuantileTable[1][m_NumberOfMatchPoints + 1] = m_ReferenceMaxValue;

  double delta = 1.0 / ( double(m_NumberOfMatchPoints) + 1.0 );
  for( unsigned int j = 1; j < m_NumberOfMatchPoints + 1; j++ )
    {
    m_QuantileTable[0][j] = m_SourceHistogram->Quantile(
        0, double(j) * delta );
    m_QuantileTable[1][j] = m_ReferenceHistogram->Quantile(
        0, double(j) * delta );
    }

  // Fill in the gradient array.
  m_Gradients.set_size( m_NumberOfMatchPoints + 1 );
  double denominator;
  for( unsigned int jj = 0; jj < m_NumberOfMatchPoints + 1; jj++ )
    {
    denominator = m_QuantileTable[0][jj + 1]
      - m_QuantileTable[0][jj];
    if( denominator != 0 )
      {
      m_Gradients[jj] = m_QuantileTable[1][jj + 1]
        - m_QuantileTable[1][jj];
      m_Gradients[jj] /= denominator;
      }
    else
      {
      m_Gradients[jj] = 0.0;
      }
    }

  denominator = m_QuantileTable[0][0] - m_SourceMinValue;
  if( denominator != 0 )
    {
    m_LowerGradient = m_QuantileTable[1][0] - m_ReferenceMinValue;
    m_LowerGradient /= denominator;
    }
  else
    {
    m_LowerGradient = 0.0;
    }

  denominator = m_QuantileTable[0][m_NumberOfMatchPoints + 1]
    - m_SourceMaxValue;
  if( denominator != 0 )
    {
    m_UpperGradient = m_QuantileTable[1][m_NumberOfMatchPoints + 1]
      - m_ReferenceMaxValue;
    m_UpperGradient /= denominator;
    }
  else
    {
    m_UpperGradient = 0.0;
    }
}

template <class TInputMesh, class TOutputMesh, class THistogramMeasurement>
void
HistogramMatchingQuadEdgeMeshFilter<TInputMesh, TOutputMesh, THistogramMeasurement>
::Transform()
{
  // Get the input and output pointers;
  InputMeshConstPointer input  = this->GetInput();
  OutputMeshPointer     output = this->GetOutput();

  InputPointDataContainerConstPointer inputDataPoint = input->GetPointData();
  OutputPointDataContainerPointer     outputDataPoint = output->GetPointData();

  // Transform the source mesh and write to output.
  typedef typename InputPointDataContainer::ConstIterator InputPointDataIterator;
  InputPointDataIterator in_Itr = inputDataPoint->Begin();

  typedef typename OutputPointDataContainer::Iterator OutputPointDataIterator;
  OutputPointDataIterator out_Itr = outputDataPoint->Begin();
  OutputPointDataIterator out_End = outputDataPoint->End();

  double srcValue, mappedValue;

  while( out_Itr != out_End )
    {
    srcValue = static_cast<double>( in_Itr.Value() );

    unsigned int j = 0;

    while( j < m_NumberOfMatchPoints + 2 )
      {
      if( srcValue < m_QuantileTable[0][j] )
        {
        break;
        }
      j++;
      }

    if( j == 0 )
      {
      // Linear interpolate from min to point[0]
      mappedValue = m_ReferenceMinValue
        + ( srcValue - m_SourceMinValue ) * m_LowerGradient;
      }
    else if( j == m_NumberOfMatchPoints + 2 )
      {
      // Linear interpolate from point[m_NumberOfMatchPoints+1] to max
      mappedValue = m_ReferenceMaxValue
        + ( srcValue - m_SourceMaxValue ) * m_UpperGradient;
      }
    else
      {
      // Linear interpolate from point[j] and point[j+1].
      mappedValue = m_QuantileTable[1][j - 1]
        + ( srcValue - m_QuantileTable[0][j - 1] ) * m_Gradients[j - 1];
      }

    out_Itr.Value() = static_cast<OutputPixelType>( mappedValue );

    ++in_Itr;
    ++out_Itr;
    }
}

template <class TInputMesh, class TOutputMesh, class THistogramMeasurement>
void
HistogramMatchingQuadEdgeMeshFilter<TInputMesh, TOutputMesh, THistogramMeasurement>
::GenerateData()
{
  this->BeforeTransform();
  this->Transform();
}

template <class TInputMesh, class TOutputMesh, class THistogramMeasurement>
void
HistogramMatchingQuadEdgeMeshFilter<TInputMesh, TOutputMesh, THistogramMeasurement>
::ComputeMinMax(
  const InputMeshType * mesh,
  THistogramMeasurement& minValue,
  THistogramMeasurement& maxValue )
{
  InputPointDataContainerConstPointer pointData = mesh->GetPointData();

  typedef typename InputPointDataContainer::ConstIterator PointDataIterator;

  PointDataIterator itr = pointData->Begin();
  PointDataIterator end = pointData->End();

  minValue = static_cast<THistogramMeasurement>( itr.Value() );
  maxValue = minValue;

  while( itr != end )
    {
    const THistogramMeasurement value = static_cast<THistogramMeasurement>( itr.Value() );

    if( value < minValue )
      {
      minValue = value;
      }
    if( value > maxValue )
      {
      maxValue = value;
      }

    ++itr;
    }
}

template <class TInputMesh, class TOutputMesh, class THistogramMeasurement>
void
HistogramMatchingQuadEdgeMeshFilter<TInputMesh, TOutputMesh, THistogramMeasurement>
::ConstructHistogram(
  const InputMeshType * mesh,
  HistogramType  * histogram,
  const THistogramMeasurement minValue,
  const THistogramMeasurement maxValue )
{
    {
    // allocate memory for the histogram
    typename HistogramType::SizeType size;
    typename HistogramType::MeasurementVectorType lowerBound;
    typename HistogramType::MeasurementVectorType upperBound;

    size.SetSize(1);
    lowerBound.SetSize(1);
    upperBound.SetSize(1);
    histogram->SetMeasurementVectorSize(1);

    size[0] = m_NumberOfHistogramLevels;
    lowerBound.Fill(minValue);
    upperBound.Fill(maxValue);

    // Initialize with equally spaced bins.
    histogram->Initialize( size, lowerBound, upperBound );
    histogram->SetToZero();
    }

  typename HistogramType::MeasurementVectorType measurement;

  measurement.SetSize(1);

  typedef typename HistogramType::MeasurementType MeasurementType;
  measurement[0] = NumericTraits<MeasurementType>::Zero;

    {
    // put each mesh scalar into the histogram
    typedef typename InputPointDataContainer::ConstIterator PointDataIterator;
    PointDataIterator itr = mesh->GetPointData()->Begin();
    PointDataIterator end = mesh->GetPointData()->End();

    while( itr != end )
      {
      InputPixelType value = itr.Value();

      if( static_cast<double>(value) >= minValue &&
          static_cast<double>(value) <= maxValue )
        {
        // add sample to histogram
        measurement[0] = value;
        histogram->IncreaseFrequency( measurement, 1 );
        }
      ++itr;
      }
    }
}

template <class TInputMesh, class TOutputMesh, class THistogramMeasurement>
void
HistogramMatchingQuadEdgeMeshFilter<TInputMesh, TOutputMesh, THistogramMeasurement>
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);

  os << indent << "NumberOfHistogramLevels: ";
  os << m_NumberOfHistogramLevels << std::endl;
  os << indent << "NumberOfMatchPoints: ";
  os << m_NumberOfMatchPoints << std::endl;

  os << indent << "Source histogram: ";
  os << m_SourceHistogram.GetPointer() << std::endl;
  os << indent << "Reference histogram: ";
  os << m_ReferenceHistogram.GetPointer() << std::endl;
  os << indent << "Output histogram: ";
  os << m_OutputHistogram.GetPointer() << std::endl;
  os << indent << "QuantileTable: " << std::endl;
  os << m_QuantileTable << std::endl;
  os << indent << "Gradients: " << std::endl;
  os << m_Gradients << std::endl;
  os << indent << "LowerGradient: ";
  os << m_LowerGradient << std::endl;
  os << indent << "UpperGradient: ";
  os << m_UpperGradient << std::endl;
}
} // end namespace itk

#endif
