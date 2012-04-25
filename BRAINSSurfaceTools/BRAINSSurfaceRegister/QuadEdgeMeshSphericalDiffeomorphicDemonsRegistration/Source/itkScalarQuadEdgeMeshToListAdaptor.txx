/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkScalarQuadEdgeMeshToListAdaptor.h,v $
  Language:  C++
  Date:      $Date: 2008-05-14 09:26:05 $
  Version:   $Revision: 1.21 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkScalarQuadEdgeMeshToListAdaptor_txx
#define __itkScalarQuadEdgeMeshToListAdaptor_txx

#include "itkScalarQuadEdgeMeshToListAdaptor.h"
#include "itkProgressReporter.h"

namespace itk
{
template <class TInputMesh>
ScalarQuadEdgeMeshToListAdaptor<TInputMesh>
::ScalarQuadEdgeMeshToListAdaptor()
{
  this->SetNumberOfRequiredInputs( 1 );
  this->SetNumberOfRequiredOutputs( 0 );
  this->SetNumberOfOutputs( 0 );

  typename ListSampleType::Pointer sample = ListSampleType::New();

  m_Sample = static_cast<ListSampleType *>( sample.GetPointer() );

  m_Sample->SetMeasurementVectorSize( 1 );
}

template <class TInputMesh>
ScalarQuadEdgeMeshToListAdaptor<TInputMesh>
::~ScalarQuadEdgeMeshToListAdaptor()
{
}

template <class TInputMesh>
void
ScalarQuadEdgeMeshToListAdaptor<TInputMesh>
::SetMesh( const InputMeshType * mesh )
{
  itkDebugMacro("setting input mesh to " << mesh);
  if( mesh != static_cast<const InputMeshType *>(this->ProcessObject::GetInput( 0 ) ) )
    {
    this->ProcessObject::SetNthInput(0, const_cast<InputMeshType *>( mesh ) );
    this->Modified();
    }
}

template <class TInputMesh>
const typename
ScalarQuadEdgeMeshToListAdaptor<TInputMesh>::InputMeshType
* ScalarQuadEdgeMeshToListAdaptor<TInputMesh>
::GetMesh() const
  {
  Self *                surrogate = const_cast<Self *>( this );
  const InputMeshType * inputMesh =
    static_cast<const InputMeshType *>( surrogate->ProcessObject::GetInput(0) );
  return inputMesh;
  }

template <class TInputMesh>
void
ScalarQuadEdgeMeshToListAdaptor<TInputMesh>
::Compute()
{
  const InputMeshType * inputMesh = this->GetMesh();

  const InputPointDataContainer * pointDataContainer = inputMesh->GetPointData();

  const unsigned int numberOfPoints = inputMesh->GetNumberOfPoints();

  ProgressReporter progress(this, 0, numberOfPoints);

  typedef typename InputPointDataContainer::ConstIterator PointDataConstIterator;

  PointDataConstIterator pointDataItr = pointDataContainer->Begin();
  PointDataConstIterator pointDataEnd = pointDataContainer->End();

  MeasurementVectorType mv;

  while( pointDataItr != pointDataEnd )
    {
    mv[0] = ( MeasurementType ) pointDataItr.Value();
    m_Sample->PushBack(mv);

    ++pointDataItr;
    }
}
} // end namespace itk

#endif
