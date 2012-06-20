/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkRescaleScalarsQuadEdgeMeshFilter.h,v $
  Language:  C++
  Date:      $Date: 2010-05-14 09:26:05 $
  Version:   $Revision: 1.21 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkRescaleScalarsQuadEdgeMeshFilter_txx
#define __itkRescaleScalarsQuadEdgeMeshFilter_txx

#include "itkRescaleScalarsQuadEdgeMeshFilter.h"
#include "itkProgressReporter.h"
#include "itkNumericTraits.h"

namespace itk
{
template <class TMesh>
RescaleScalarsQuadEdgeMeshFilter<TMesh>
::RescaleScalarsQuadEdgeMeshFilter()
{
  this->SetNumberOfRequiredInputs( 1 );
  this->SetNumberOfRequiredOutputs( 1 );
  this->SetNumberOfIndexedOutputs( 1 );

  this->SetNthOutput( 0, OutputMeshType::New() );

  this->m_OutputMaximum = NumericTraits<OutputPixelType>::max();
  this->m_OutputMinimum = NumericTraits<OutputPixelType>::NonpositiveMin();

  this->m_InputMaximum = NumericTraits<InputPixelType>::Zero;
  this->m_InputMinimum = NumericTraits<InputPixelType>::max();
}

template <class TMesh>
RescaleScalarsQuadEdgeMeshFilter<TMesh>
::~RescaleScalarsQuadEdgeMeshFilter()
{
}

template <class TMesh>
void
RescaleScalarsQuadEdgeMeshFilter<TMesh>
::SetInputMesh( const InputMeshType * mesh )
{
  itkDebugMacro("setting input mesh to " << mesh);
  if( mesh != static_cast<const InputMeshType *>(this->ProcessObject::GetInput( 0 ) ) )
    {
    this->ProcessObject::SetNthInput(0, const_cast<InputMeshType *>( mesh ) );
    this->Modified();
    }
}

template <class TMesh>
const typename
RescaleScalarsQuadEdgeMeshFilter<TMesh>::InputMeshType
* RescaleScalarsQuadEdgeMeshFilter<TMesh>
::GetInputMesh() const
  {
  Self *                surrogate = const_cast<Self *>( this );
  const InputMeshType * inputMesh =
    static_cast<const InputMeshType *>( surrogate->ProcessObject::GetInput(0) );

  return inputMesh;
  }

template <class TMesh>
void
RescaleScalarsQuadEdgeMeshFilter<TMesh>
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

  if( this->m_InputMinimum != this->m_InputMaximum )
    {
    this->m_Scale =
      (static_cast<double>( this->m_OutputMaximum )
       - static_cast<double>( this->m_OutputMinimum ) )
      / (static_cast<double>( this->m_InputMaximum )
         - static_cast<double>( this->m_InputMinimum ) );
    }
  else if( this->m_InputMaximum != NumericTraits<InputPixelType>::Zero )
    {
    this->m_Scale =
      (static_cast<double>( this->m_OutputMaximum )
       - static_cast<double>( this->m_OutputMinimum ) )
      / static_cast<double>( this->m_InputMaximum );
    }
  else
    {
    this->m_Scale = 0.0;
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
    pointDataItr.Value() = static_cast<OutputPixelType>(
        this->m_OutputMinimum + ( (pointDataItr.Value() - this->m_InputMinimum)
                                  * this->m_Scale) );

    ++pointDataItr;
    }
}
} // end namespace itk

#endif
