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
#ifndef __itkQuadEdgeMeshSplitFilter_hxx
#define __itkQuadEdgeMeshSplitFilter_hxx

#include "itkQuadEdgeMeshSplitFilter.h"

namespace itk
{
template <typename TInputMesh, typename TOutputMesh>
QuadEdgeMeshSplitFilter<TInputMesh, TOutputMesh>::QuadEdgeMeshSplitFilter() : Superclass(), m_StartCellId( 0 ), m_Area(
    0. )
{
  this->Superclass::SetNumberOfRequiredInputs( 1 );
  this->SetNumberOfRequiredOutputs( 2 );
  this->SetNumberOfIndexedOutputs( 2 );
  this->SetNthOutput( 0, OutputMeshType::New() );
  this->SetNthOutput( 1, OutputMeshType::New() );
}

template <typename TInputMesh, typename TOutputMesh>
QuadEdgeMeshSplitFilter<TInputMesh, TOutputMesh>::
~QuadEdgeMeshSplitFilter()
{
}

template <typename TInputMesh, typename TOutputMesh>
void
QuadEdgeMeshSplitFilter<TInputMesh, TOutputMesh>::GenerateData()
{
  InputMeshConstPointer input = this->GetInput();

  OutputMeshPointer output0 = this->GetOutput( 0 );

  output0->SetCellsAllocationMethod(
    OutputMeshType::CellsAllocatedDynamicallyCellByCell );
  OutputMeshPointer output1 = this->GetOutput( 1 );
  output1->SetCellsAllocationMethod(
    OutputMeshType::CellsAllocatedDynamicallyCellByCell );

  InputCellsContainerConstPointer cells = input->GetCells();

  FMMSeedVectorType seeds;

  seeds.push_back( m_SeedFaces[0] );
  seeds.push_back( m_SeedFaces[1] );

  FMMPointer fast_marching_clustering = FMMType::New();
  fast_marching_clustering->SetMesh( input );
  fast_marching_clustering->SetSeedFaces( seeds );

  InputPolygonType*  poly;
  FMMClusterType     cluster = fast_marching_clustering->Evaluate( 0 );
  FMMClusterIterator cluster_it;

  std::list<OutputPointIdList> list_face;
  for( cluster_it = cluster.begin();
       cluster_it != cluster.end();
       ++cluster_it )
    {
    poly = dynamic_cast<InputPolygonType *>(
        cells->GetElement( *cluster_it ) );
    list_face.push_back( AddFacePointsToOutputMesh( output0, poly ) );
    }

  typename std::list<OutputPointIdList>::iterator list_face_it;
  for( list_face_it = list_face.begin();
       list_face_it != list_face.end();
       ++list_face_it )
    {
    output0->AddFaceWithSecurePointList( *list_face_it );
    }
  list_face.clear();

  cluster = fast_marching_clustering->Evaluate( 1 );
  for( cluster_it = cluster.begin();
       cluster_it != cluster.end();
       ++cluster_it )
    {
    poly = dynamic_cast<InputPolygonType *>(
        cells->GetElement( *cluster_it ) );
    list_face.push_back( AddFacePointsToOutputMesh( output1, poly ) );
    }
  for( list_face_it = list_face.begin();
       list_face_it != list_face.end();
       ++list_face_it )
    {
    output1->AddFaceWithSecurePointList( *list_face_it );
    }
}

template <typename TInputMesh, typename TOutputMesh>
typename QuadEdgeMeshSplitFilter<TInputMesh, TOutputMesh>::OutputPointIdList
QuadEdgeMeshSplitFilter<TInputMesh, TOutputMesh>::AddFacePointsToOutputMesh( OutputMeshType* ioMesh,
                                                                             InputPolygonType* iPoly )
{
  InputMeshConstPointer input = this->GetInput();

  InputQEType*         edge = iPoly->GetEdgeRingEntry();
  InputQEType*         temp = edge;
  InputPointType       p;
  OutputPointType      q;
  InputPointIdentifier id_org;
  OutputPointIdList    oPoints;

  OutputPointsContainerPointer points_container = ioMesh->GetPoints();

  do
    {
    id_org = temp->GetOrigin();
    oPoints.push_back( id_org );
    if( !points_container->GetElementIfIndexExists( id_org, &q ) )
      {
      p = input->GetPoint( id_org );
      q.CastFrom( p );
      ioMesh->SetPoint( id_org, q );
      }
    temp = temp->GetLnext();
    }
  while( temp != edge );

  return oPoints;
}
}

#endif
