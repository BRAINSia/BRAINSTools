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
#ifndef __itkHistogramMatchingQuadEdgeMeshFilter_hxx
#define __itkHistogramMatchingQuadEdgeMeshFilter_hxx

#include "itkHistogramMatchingQuadEdgeMeshFilter.h"
#include "itkNumericTraits.h"
#include <vector>

namespace itk
{
template <typename TInputMesh, typename TOutputMesh, typename THistogramMeasurement>
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

template <typename TInputMesh, typename TOutputMesh, typename THistogramMeasurement>
HistogramMatchingQuadEdgeMeshFilter<TInputMesh, TOutputMesh, THistogramMeasurement>
::~HistogramMatchingQuadEdgeMeshFilter()
{
}

template <typename TInputMesh, typename TOutputMesh, typename THistogramMeasurement>
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

template <typename TInputMesh, typename TOutputMesh, typename THistogramMeasurement>
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

template <typename TInputMesh, typename TOutputMesh, typename THistogramMeasurement>
void
HistogramMatchingQuadEdgeMeshFilter<TInputMesh, TOutputMesh, THistogramMeasurement>
::SetSourceMesh( const InputMeshType * source )
{
  this->SetInput( source );
}

template <typename TInputMesh, typename TOutputMesh, typename THistogramMeasurement>
const typename HistogramMatchingQuadEdgeMeshFilter<TInputMesh, TOutputMesh, THistogramMeasurement>::InputMeshType
* HistogramMatchingQuadEdgeMeshFilter<TInputMesh, TOutputMesh, THistogramMeasurement>
::GetSourceMesh( void ) const
  {
  return this->GetInput();
  }

template <typename TInputMesh, typename TOutputMesh, typename THistogramMeasurement>
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

template <typename TInputMesh, typename TOutputMesh, typename THistogramMeasurement>
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
  using InputPointDataIterator = typename InputPointDataContainer::ConstIterator;
  InputPointDataIterator in_Itr = inputDataPoint->Begin();

  using OutputPointDataIterator = typename OutputPointDataContainer::Iterator;
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

template <typename TInputMesh, typename TOutputMesh, typename THistogramMeasurement>
void
HistogramMatchingQuadEdgeMeshFilter<TInputMesh, TOutputMesh, THistogramMeasurement>
::GenerateData()
{
  this->BeforeTransform();
  this->Transform();
}

template <typename TInputMesh, typename TOutputMesh, typename THistogramMeasurement>
void
HistogramMatchingQuadEdgeMeshFilter<TInputMesh, TOutputMesh, THistogramMeasurement>
::ComputeMinMax(
  const InputMeshType * mesh,
  THistogramMeasurement& minValue,
  THistogramMeasurement& maxValue )
{
  InputPointDataContainerConstPointer pointData = mesh->GetPointData();

  using PointDataIterator = typename InputPointDataContainer::ConstIterator;

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

template <typename TInputMesh, typename TOutputMesh, typename THistogramMeasurement>
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

  using MeasurementType = typename HistogramType::MeasurementType;
  measurement[0] = NumericTraits<MeasurementType>::ZeroValue();

    {
    // put each mesh scalar into the histogram
    using PointDataIterator = typename InputPointDataContainer::ConstIterator;
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
        histogram->IncreaseFrequencyOfMeasurement( measurement, 1 );
        }
      ++itr;
      }
    }
}

template <typename TInputMesh, typename TOutputMesh, typename THistogramMeasurement>
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
