/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkNodeScalarGradientCalculator.h,v $
  Language:  C++
  Date:      $Date: 2008-05-14 09:26:05 $
  Version:   $Revision: 1.21 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkAssignScalarValuesQuadEdgeMeshFilter_txx
#define __itkAssignScalarValuesQuadEdgeMeshFilter_txx

#include "itkAssignScalarValuesQuadEdgeMeshFilter.h"
#include "itkProgressReporter.h"
#include "itkNumericTraitsVectorPixel.h"

namespace itk
{
template <class TInputMesh, class TSourceMesh, class TOutputMesh>
AssignScalarValuesQuadEdgeMeshFilter<TInputMesh, TSourceMesh, TOutputMesh>
::AssignScalarValuesQuadEdgeMeshFilter()
{
  this->SetNumberOfRequiredInputs( 2 );
  this->SetNumberOfRequiredOutputs( 1 );
  this->SetNumberOfIndexedOutputs( 1 );

  this->SetNthOutput( 0, OutputMeshType::New() );
}

template <class TInputMesh, class TSourceMesh, class TOutputMesh>
AssignScalarValuesQuadEdgeMeshFilter<TInputMesh, TSourceMesh, TOutputMesh>
::~AssignScalarValuesQuadEdgeMeshFilter()
{
}

template <class TInputMesh, class TSourceMesh, class TOutputMesh>
void
AssignScalarValuesQuadEdgeMeshFilter<TInputMesh, TSourceMesh, TOutputMesh>
::SetSourceMesh( const SourceMeshType * mesh )
{
  itkDebugMacro("setting source Mesh to " << mesh);
  if( mesh != static_cast<const SourceMeshType *>(this->ProcessObject::GetInput( 0 ) ) )
    {
    this->ProcessObject::SetNthInput(1, const_cast<SourceMeshType *>( mesh ) );
    this->Modified();
    }
}

template <class TInputMesh, class TSourceMesh, class TOutputMesh>
const typename
AssignScalarValuesQuadEdgeMeshFilter<TInputMesh, TSourceMesh, TOutputMesh>::SourceMeshType
* AssignScalarValuesQuadEdgeMeshFilter<TInputMesh, TSourceMesh, TOutputMesh>
::GetSourceMesh() const
  {
  Self *                 surrogate = const_cast<Self *>( this );
  const SourceMeshType * sourceMesh =
    static_cast<const SourceMeshType *>( surrogate->ProcessObject::GetInput(1) );
  return sourceMesh;
  }

template <class TInputMesh, class TSourceMesh, class TOutputMesh>
void
AssignScalarValuesQuadEdgeMeshFilter<TInputMesh, TSourceMesh, TOutputMesh>
::SetInputMesh( const InputMeshType * mesh )
{
  itkDebugMacro("setting input Mesh to " << mesh);
  if( mesh != static_cast<const InputMeshType *>(this->ProcessObject::GetInput( 1 ) ) )
    {
    this->ProcessObject::SetNthInput(0, const_cast<InputMeshType *>( mesh ) );
    this->Modified();
    }
}

template <class TInputMesh, class TSourceMesh, class TOutputMesh>
const typename
AssignScalarValuesQuadEdgeMeshFilter<TInputMesh, TSourceMesh, TOutputMesh>::InputMeshType
* AssignScalarValuesQuadEdgeMeshFilter<TInputMesh, TSourceMesh, TOutputMesh>
::GetInputMesh() const
  {
  Self *                surrogate = const_cast<Self *>( this );
  const InputMeshType * inputMesh =
    static_cast<const InputMeshType *>( surrogate->ProcessObject::GetInput(0) );
  return inputMesh;
  }

template <class TInputMesh, class TSourceMesh, class TOutputMesh>
void
AssignScalarValuesQuadEdgeMeshFilter<TInputMesh, TSourceMesh, TOutputMesh>
::GenerateData()
{
  this->CopyInputMeshToOutputMesh();

  const SourceMeshType * sourceMesh = this->GetSourceMesh();

  const SourcePointDataContainer * sourcePointData = sourceMesh->GetPointData();

  OutputPointDataContainerPointer outputPointData = this->GetOutput()->GetPointData();

  if( outputPointData.IsNull() )
    {
    outputPointData = OutputPointDataContainer::New();
    this->GetOutput()->SetPointData( outputPointData );
    }

  const unsigned int numberOfNodes = sourceMesh->GetNumberOfPoints();

  outputPointData->Reserve( numberOfNodes );

  if( !sourcePointData )
    {
    itkExceptionMacro("Source PointData is missing");
    }

  OutputPointDataContainerIterator outputDataItr = outputPointData->Begin();

  typedef typename SourcePointDataContainer::ConstIterator SourcePointDataIterator;
  SourcePointDataIterator sourceDataItr = sourcePointData->Begin();
  SourcePointDataIterator sourceDataEnd = sourcePointData->End();

  while( sourceDataItr != sourceDataEnd )
    {
    outputDataItr.Value() = sourceDataItr.Value();

    ++outputDataItr;
    ++sourceDataItr;
    }
}
} // end namespace itk

#endif
