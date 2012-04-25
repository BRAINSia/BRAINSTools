#ifndef __itkQuadEdgeMeshSplitFilter_txx
#define __itkQuadEdgeMeshSplitFilter_txx

#include "itkQuadEdgeMeshSplitFilter.txx"

namespace itk
{
template <class TInputMesh, class TOutputMesh>
QuadEdgeMeshSplitFilter<TInputMesh, TOutputMesh>::QuadEdgeMeshSplitFilter() : Superclass(), m_StartCellId( 0 ), m_Area(
    0. )
{
  this->Superclass::SetNumberOfRequiredInputs( 1 );
  this->SetNumberOfRequiredOutputs( 2 );
  this->SetNumberOfOutputs( 2 );
  this->SetNthOutput( 0, OutputMeshType::New() );
  this->SetNthOutput( 1, OutputMeshType::New() );
}

template <class TInputMesh, class TOutputMesh>
QuadEdgeMeshSplitFilter<TInputMesh, TOutputMesh>::
~QuadEdgeMeshSplitFilter()
{
}

template <class TInputMesh, class TOutputMesh>
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

template <class TInputMesh, class TOutputMesh>
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
