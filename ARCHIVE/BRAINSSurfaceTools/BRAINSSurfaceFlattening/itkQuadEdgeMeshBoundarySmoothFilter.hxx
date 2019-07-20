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
#ifndef __itkQuadEdgeMeshBoundarySmoothFilter_hxx
#define __itkQuadEdgeMeshBoundarySmoothFilter_hxx

#include "itkQuadEdgeMeshBoundarySmoothFilter.h"
#include "itkQuadEdgeMeshBoundaryEdgesMeshFunction.h"

namespace itk
{
template < typename TInputMesh, typename TOutputMesh >
QuadEdgeMeshBoundarySmoothFilter< TInputMesh, TOutputMesh >::QuadEdgeMeshBoundarySmoothFilter()
{
  this->SetNumberOfRequiredInputs( 2 );
  this->SetNumberOfRequiredOutputs( 2 );
  this->SetNumberOfIndexedOutputs( 2 );
  this->SetNthOutput( 0, OutputMeshType::New() );
  this->SetNthOutput( 1, OutputMeshType::New() );

  m_Iterations = 10;
}

template < typename TInputMesh, typename TOutputMesh >
QuadEdgeMeshBoundarySmoothFilter< TInputMesh, TOutputMesh >::~QuadEdgeMeshBoundarySmoothFilter()
{}

template < typename TInputMesh, typename TOutputMesh >
void
QuadEdgeMeshBoundarySmoothFilter< TInputMesh, TOutputMesh >::SetInputMesh1( const InputMeshType * mesh1 )
{
  itkDebugMacro( "setting input mesh 1 to " << mesh1 );
  if ( mesh1 != static_cast< const InputMeshType * >( this->ProcessObject::GetInput( 0 ) ) )
  {
    this->ProcessObject::SetNthInput( 0, const_cast< InputMeshType * >( mesh1 ) );
    this->Modified();
  }
}

template < typename TInputMesh, typename TOutputMesh >
const typename QuadEdgeMeshBoundarySmoothFilter< TInputMesh, TOutputMesh >::InputMeshType *
QuadEdgeMeshBoundarySmoothFilter< TInputMesh, TOutputMesh >::GetInputMesh1() const
{
  Self *                surrogate = const_cast< Self * >( this );
  const InputMeshType * inputMesh1 = static_cast< const InputMeshType * >( surrogate->ProcessObject::GetInput( 0 ) );
  return inputMesh1;
}

template < typename TInputMesh, typename TOutputMesh >
void
QuadEdgeMeshBoundarySmoothFilter< TInputMesh, TOutputMesh >::SetInputMesh2( const InputMeshType * mesh2 )
{
  itkDebugMacro( "setting input mesh 2 to " << mesh2 );
  if ( mesh2 != static_cast< const InputMeshType * >( this->ProcessObject::GetInput( 1 ) ) )
  {
    this->ProcessObject::SetNthInput( 1, const_cast< InputMeshType * >( mesh2 ) );
    this->Modified();
  }
}

template < typename TInputMesh, typename TOutputMesh >
const typename QuadEdgeMeshBoundarySmoothFilter< TInputMesh, TOutputMesh >::InputMeshType *
QuadEdgeMeshBoundarySmoothFilter< TInputMesh, TOutputMesh >::GetInputMesh2() const
{
  Self *                surrogate = const_cast< Self * >( this );
  const InputMeshType * inputMesh2 = static_cast< const InputMeshType * >( surrogate->ProcessObject::GetInput( 1 ) );
  return inputMesh2;
}

template < typename TInputMesh, typename TOutputMesh >
typename QuadEdgeMeshBoundarySmoothFilter< TInputMesh, TOutputMesh >::OutputMeshType *
QuadEdgeMeshBoundarySmoothFilter< TInputMesh, TOutputMesh >::GetOutputMesh1()
{
  Self *           surrogate = const_cast< Self * >( this );
  OutputMeshType * outputMesh1 = static_cast< OutputMeshType * >( surrogate->ProcessObject::GetOutput( 0 ) );

  return outputMesh1;
}

template < typename TInputMesh, typename TOutputMesh >
typename QuadEdgeMeshBoundarySmoothFilter< TInputMesh, TOutputMesh >::OutputMeshType *
QuadEdgeMeshBoundarySmoothFilter< TInputMesh, TOutputMesh >::GetOutputMesh2()
{
  Self *           surrogate = const_cast< Self * >( this );
  OutputMeshType * outputMesh2 = static_cast< OutputMeshType * >( surrogate->ProcessObject::GetOutput( 1 ) );

  return outputMesh2;
}

template < typename TInputMesh, typename TOutputMesh >
void
QuadEdgeMeshBoundarySmoothFilter< TInputMesh, TOutputMesh >::GenerateData()
{
  this->CopyInputMeshesToOutputMeshes();
  unsigned long NumberOfCellsInput, NumberOfCellstmp;

  NumberOfCellsInput = this->GetInputMesh1()->GetNumberOfCells() + this->GetInputMesh2()->GetNumberOfCells();

  // assign which mesh is going to delete face
  OutputMeshPointer mesh1 = this->GetOutputMesh1();
  OutputMeshPointer mesh2 = this->GetOutputMesh2();
  for ( int Iter_i = 0; Iter_i < m_Iterations; Iter_i++ )
  {
    int num1 = this->AdjustBoundary( mesh1, mesh2 );
    NumberOfCellstmp = mesh1->GetNumberOfCells() + mesh2->GetNumberOfCells();
    if ( NumberOfCellstmp != NumberOfCellsInput )
    {
      itkExceptionMacro( "Number of cells is changed." );
    }

    int num2 = this->AdjustBoundary( mesh2, mesh1 );
    NumberOfCellstmp = mesh1->GetNumberOfCells() + mesh2->GetNumberOfCells();
    if ( NumberOfCellstmp != NumberOfCellsInput )
    {
      itkExceptionMacro( "Number of cells is changed." );
    }

    if ( ( num1 == 0 ) && ( num2 == 0 ) )
    {
      break;
    }
  }

  this->SetNthOutput( 0, mesh1 );
  this->SetNthOutput( 1, mesh2 );
}

template < typename TInputMesh, typename TOutputMesh >
int
QuadEdgeMeshBoundarySmoothFilter< TInputMesh, TOutputMesh >::AdjustBoundary( OutputMeshType * deleteMesh,
                                                                             OutputMeshType * addMesh )
{
  using BoundaryFunctionType = typename itk::QuadEdgeMeshBoundaryEdgesMeshFunction< TOutputMesh >;

  typename BoundaryFunctionType::Pointer boundaryFunction = BoundaryFunctionType::New();

  int numChange = 0;

  // get the boundary of deleteMesh
  EdgeListPointerType boundaryList;

  boundaryList = boundaryFunction->Evaluate( *deleteMesh );
  OutputQEPrimal * bdryEdge = boundaryList->front();
  OutputQEIterator it = bdryEdge->BeginGeomLnext();

  // go through boundary of deleteMesh
  OutputPointIdentifier         pId0, pId1, pId2;
  OutputPointIdentifier         pId0_2, pId1_2, pId2_2;
  OutputPointsContainerIterator pIt, pIt_end;

  std::vector< InputCellIdentifier > teethList;

  while ( it != bdryEdge->EndGeomLnext() )
  {
    pId0 = it.Value()->GetOrigin();
    pId1 = it.Value()->GetDestination();

    OutputQEPrimal * tmpLeft = it.Value()->GetLnext();

    pId2 = tmpLeft->GetDestination();

    OutputQEPrimal * e = deleteMesh->FindEdge( pId0, pId2 );

    // find two edges on the boundary belong to the same triangle
    if ( e != (OutputQEPrimal *)nullptr )
    {
      // add first
      pIt = addMesh->GetPoints()->Begin();
      pIt_end = addMesh->GetPoints()->End();
      int numPoints = 0;
      while ( pIt != pIt_end )
      {
        if ( pIt.Value() == deleteMesh->GetPoint( pId0 ) )
        {
          pId0_2 = pIt.Index();
          numPoints += 1;
        }
        if ( pIt.Value() == deleteMesh->GetPoint( pId1 ) )
        {
          pId1_2 = pIt.Index();
          numPoints += 1;
        }
        if ( pIt.Value() == deleteMesh->GetPoint( pId2 ) )
        {
          pId2_2 = pIt.Index();
          numPoints += 1;
        }

        if ( numPoints == 3 )
        {
          break;
        }
        pIt++;
      }

      // actual add
      OutputQEPrimal * tmpEdge = addMesh->AddFaceTriangle( pId0_2, pId1_2, pId2_2 );
      if ( tmpEdge == nullptr )
      {
        tmpEdge = addMesh->AddFaceTriangle( pId0_2, pId2_2, pId1_2 );
      }

      if ( e->IsLeftSet() )
      {
        // remember the cellId
        teethList.push_back( e->GetLeft() );
      }

      numChange += 1;
    }
    it++;
  }

  boundaryList->clear();

  // create a new mesh
  // add the cells from deleteMesh to newMesh,
  // avoid the cells in teethList.

  OutputMeshPointer newMesh = OutputMeshType::New();
  newMesh->SetCellsAllocationMethod( OutputMeshType::CellsAllocatedDynamicallyCellByCell );

  OutputCellsContainerConstPointer cells = deleteMesh->GetCells();

  OutputCellsContainerConstIterator cells_it = cells->Begin();
  OutputCellsContainerConstIterator cells_end = cells->End();

  OutputPolygonCellType * poly;
  OutputPointIdListType   securePointIdList;

  while ( cells_it != cells_end )
  {
    if ( std::find( teethList.begin(), teethList.end(), cells_it.Index() ) == teethList.end() )
    {
      // not a teeth
      poly = dynamic_cast< InputPolygonCellType * >( cells->GetElement( cells_it.Index() ) );

      // add points first
      OutputQEPrimal *      firstEdge = poly->GetEdgeRingEntry();
      OutputQEPrimal *      temp = firstEdge;
      OutputPointType       p;
      OutputPointType       q;
      OutputPointIdentifier id_org;

      OutputPointsContainerPointer points_container = newMesh->GetPoints();

      do
      {
        id_org = temp->GetOrigin();
        securePointIdList.push_back( id_org );
        if ( !points_container->GetElementIfIndexExists( id_org, &p ) )
        {
          p = deleteMesh->GetPoint( id_org );
          q.CastFrom( p );
          newMesh->SetPoint( id_org, q );
        }
        temp = temp->GetLnext();
      } while ( temp != firstEdge );

      newMesh->AddFaceWithSecurePointList( securePointIdList );
      securePointIdList.clear();
    }

    cells_it++;
  }

  deleteMesh->Clear();
  CopyMeshToMesh< OutputMeshType, OutputMeshType >( newMesh, deleteMesh );

  return numChange;
}

template < typename TInputMesh, typename TOutputMesh >
void
QuadEdgeMeshBoundarySmoothFilter< TInputMesh, TOutputMesh >::CopyInputMeshesToOutputMeshes()
{
  const InputMeshType * in1 = this->GetInputMesh1();
  OutputMeshType *      out1 = this->GetOutputMesh1();

  CopyMeshToMesh( in1, out1 );

  const InputMeshType * in2 = this->GetInputMesh2();
  OutputMeshType *      out2 = this->GetOutputMesh2();

  CopyMeshToMesh( in2, out2 );
}
} // namespace itk

#endif
