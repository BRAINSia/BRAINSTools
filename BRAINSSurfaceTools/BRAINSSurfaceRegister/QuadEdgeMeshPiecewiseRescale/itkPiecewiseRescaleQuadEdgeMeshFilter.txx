/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkPiecewiseRescaleQuadEdgeMeshFilter.h,v $
  Language:  C++
  Date:      $Date: 2011-05-03 13:09:05 $
  Version:   $Revision: 1.0 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkPiecewiseRescaleQuadEdgeMeshFilter_txx
#define __itkPiecewiseRescaleQuadEdgeMeshFilter_txx

#include "itkPiecewiseRescaleQuadEdgeMeshFilter.h"
#include "itkProgressReporter.h"
#include "itkNumericTraits.h"

namespace itk
{
template <class TInputMesh, class TOutputMesh>
PiecewiseRescaleQuadEdgeMeshFilter<TInputMesh, TOutputMesh>
::PiecewiseRescaleQuadEdgeMeshFilter()
{
  this->SetNumberOfRequiredInputs( 1 );
  this->SetNumberOfRequiredOutputs( 1 );
  this->SetNumberOfOutputs( 1 );

  this->SetNthOutput( 0, OutputMeshType::New() );

  this->m_OutputMaximum = NumericTraits<OutputPixelType>::max();
  this->m_OutputMinimum = NumericTraits<OutputPixelType>::NonpositiveMin();

  this->m_cValue = NumericTraits<InputPixelType>::Zero;
  this->m_InputMaximum = NumericTraits<InputPixelType>::min();
  this->m_InputMinimum = NumericTraits<InputPixelType>::max();
}

template <class TInputMesh, class TOutputMesh>
PiecewiseRescaleQuadEdgeMeshFilter<TInputMesh, TOutputMesh>
::~PiecewiseRescaleQuadEdgeMeshFilter()
{
}

template <class TInputMesh, class TOutputMesh>
void
PiecewiseRescaleQuadEdgeMeshFilter<TInputMesh, TOutputMesh>
::SetInputMesh( const InputMeshType * mesh )
{
  itkDebugMacro("setting input mesh to " << mesh);
  if( mesh != static_cast<const InputMeshType *>(this->ProcessObject::GetInput( 0 ) ) )
    {
    this->ProcessObject::SetNthInput(0, const_cast<InputMeshType *>( mesh ) );
    this->Modified();
    }
}

template <class TInputMesh, class TOutputMesh>
const typename
PiecewiseRescaleQuadEdgeMeshFilter<TInputMesh, TOutputMesh>::InputMeshType
* PiecewiseRescaleQuadEdgeMeshFilter<TInputMesh, TOutputMesh>
::GetInputMesh() const
  {
  Self *                surrogate = const_cast<Self *>( this );
  const InputMeshType * inputMesh =
    static_cast<const InputMeshType *>( surrogate->ProcessObject::GetInput(0) );

  return inputMesh;
  }

template <class TInputMesh, class TOutputMesh>
void
PiecewiseRescaleQuadEdgeMeshFilter<TInputMesh, TOutputMesh>
::GenerateData()
{
  // Copy the input mesh into the output mesh.
  this->CopyInputMeshToOutputMesh();

  OutputMeshType * outputMesh = this->GetOutput();

  const InputMeshType * inputMesh = this->GetInputMesh();

  //
  // Visit all nodes of the input mesh
  //
  InputPointDataContainerConstPointer inputPointData = inputMesh->GetPointData();
  typedef typename InputPointDataContainer::ConstIterator InputPointDataIterator;

  InputPointDataIterator inputPointDataItr = inputPointData->Begin();
  InputPointDataIterator inputPointDataEnd = inputPointData->End();

  while( inputPointDataItr != inputPointDataEnd )
    {
    if( inputPointDataItr.Value() > this->m_InputMaximum )
      {
      this->m_InputMaximum = inputPointDataItr.Value();
      }

    if( inputPointDataItr.Value() < this->m_InputMinimum )
      {
      this->m_InputMinimum = inputPointDataItr.Value();
      }

    ++inputPointDataItr;
    }

  // verify the cValue
  if( m_cValue > m_InputMaximum || m_cValue < m_InputMinimum )
    {
    itkExceptionMacro("cValue should be in [InputMinimum, InputMaximum]!");
    }

  if( (this->m_InputMinimum != this->m_cValue) &&
      (this->m_InputMaximum != this->m_cValue) )
    {
    this->m_Scale_a =
      (static_cast<double>( this->m_cValue )
       - static_cast<double>( this->m_OutputMinimum ) )
      / (static_cast<double>( this->m_cValue )
         - static_cast<double>( this->m_InputMinimum ) );

    this->m_Scale_b =
      (static_cast<double>( this->m_OutputMaximum )
       - static_cast<double>( this->m_cValue ) )
      / (static_cast<double>( this->m_InputMaximum )
         - static_cast<double>( this->m_cValue ) );
    }
  else if( this->m_InputMinimum != this->m_InputMaximum )
    {
    this->m_Scale_a =
      (static_cast<double>( this->m_OutputMaximum )
       - static_cast<double>( this->m_OutputMinimum ) )
      / (static_cast<double>( this->m_InputMaximum )
         - static_cast<double>( this->m_InputMinimum ) );
    this->m_Scale_b = this->m_Scale_a;
    }
  else if( this->m_InputMaximum != NumericTraits<InputPixelType>::Zero )
    {
    this->m_Scale_a =
      (static_cast<double>( this->m_OutputMaximum )
       - static_cast<double>( this->m_OutputMinimum ) )
      / static_cast<double>( this->m_InputMaximum );
    this->m_Scale_b = this->m_Scale_a;
    }
  else
    {
    this->m_Scale_a = 0.0;
    this->m_Scale_b = 0.0;
    }

  //
  // Visit all nodes of the output mesh
  // and change the scalars
  OutputPointDataContainer * pointData = outputMesh->GetPointData();

  typedef typename OutputPointDataContainer::Iterator PointDataIterator;

  PointDataIterator pointDataItr = pointData->Begin();
  PointDataIterator pointDataEnd = pointData->End();

  while( pointDataItr != pointDataEnd )
    {
    if( pointDataItr.Value() < static_cast<OutputPixelType>(m_cValue) )
      {
      pointDataItr.Value() = static_cast<OutputPixelType>(
          this->m_OutputMinimum + ( (pointDataItr.Value() - this->m_InputMinimum)
                                    * this->m_Scale_a) );
      }
    else if( pointDataItr.Value() > static_cast<OutputPixelType>(m_cValue) )
      {
      pointDataItr.Value() = static_cast<OutputPixelType>(
          this->m_cValue + ( (pointDataItr.Value() - this->m_cValue)
                             * this->m_Scale_b) );
      }
    else
      {
      pointDataItr.Value() = static_cast<OutputPixelType>(m_cValue);
      }

    ++pointDataItr;
    }
}
} // end namespace itk

#endif
