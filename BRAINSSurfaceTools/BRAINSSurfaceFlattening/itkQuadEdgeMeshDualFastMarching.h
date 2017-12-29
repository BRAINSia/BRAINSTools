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
#ifndef __itkQuadEdgeMeshDualFastMarching_h
#define __itkQuadEdgeMeshDualFastMarching_h

#include <itkPriorityQueueContainer.h>
#include <itkTriangleHelper.h>
#include <itkQuadEdgeMeshFunctionBase.h>
#include <itkQuadEdgeMeshPolygonCell.h>
#include <vector>
#include <list>

#include "itkQuadEdgeMeshDualMetric.h"

namespace itk
{
/**
* \class QuadEdgeMeshDualFastMarching
* \brief Fast marching implemented in a Dijkstra-like way. Given a list of
* seeds (faces) compute clusters (one face belongs to one cluster if its
* "distance" to the seed is smaller the distance to the other one).
* \note The distance between 2 faces is defined by TMetric.
* */
template <class TMesh,
          class TMetric = QuadEdgeMeshDualSquaredEuclideanMetric<TMesh> >
class QuadEdgeMeshDualFastMarching :
  public QuadEdgeMeshFunctionBase<TMesh,
                                  std::list<typename TMesh::CellIdentifier> >
{
public:
  typedef QuadEdgeMeshDualFastMarching Self;
  typedef QuadEdgeMeshFunctionBase<
      TMesh,
      std::list<typename TMesh::CellIdentifier> >
    Superclass;
  typedef SmartPointer<Self>       Pointer;
  typedef SmartPointer<const Self> ConstPointer;

  itkNewMacro( Self );

  itkTypeMacro( QuadEdgeMeshDualFastMarching, QuadEdgeMeshFunctionBase );

  typedef TMesh                              MeshType;
  typedef typename MeshType::ConstPointer    MeshConstPointer;
  typedef typename MeshType::PointType       PointType;
  typedef typename MeshType::CellType        CellType;
  typedef typename MeshType::CellIdentifier  CellIdentifier;
  typedef typename MeshType::CellAutoPointer CellAutoPointer;
  typedef typename MeshType::CellsContainerConstPointer
    CellsContainerConstPointer;
  typedef typename MeshType::CellsContainerConstIterator
    CellsContainerConstIterator;
  typedef typename MeshType::PointIdIterator PointIdIterator;
  typedef typename MeshType::PointIdentifier PointIdentifier;
  typedef typename MeshType::QEType          QEType;

  typedef TriangleHelper<PointType> TriangleType;

  typedef QuadEdgeMeshPolygonCell<CellType>     PolygonType;
  typedef typename PolygonType::SelfAutoPointer PolygonAutoPointer;

  typedef TMetric                        MetricType;
  typedef typename MetricType::ValueType MetricValueType;

  typedef typename Superclass::OutputType CellIdentifierListType;

  typedef CellIdentifierListType               ClusterType;
  typedef typename ClusterType::iterator       ClusterIterator;
  typedef std::vector<ClusterType>             ClusterVectorType;
  typedef typename ClusterVectorType::iterator ClusterVectorIterator;

  typedef std::vector<CellIdentifier>       SeedVectorType;
  typedef typename SeedVectorType::iterator SeedVectorIterator;
public:

  struct Item
    {
    Item() : m_Face( 0 ), m_Value( -1. ), m_Cluster( 0 )
    {
    }

    Item( const CellIdentifier& iFace,
          const MetricValueType& iValue,
          const CellIdentifier& iCluster ) :
      m_Face( iFace ), m_Value( iValue ), m_Cluster( iCluster )
    {
    }

    Item(int x) : m_Face( 0 ), m_Value( -1. ), m_Cluster( x )
    {
    }

    ~Item()
    {
    }

    CellIdentifier m_Face;
    MetricValueType m_Value;
    CellIdentifier m_Cluster;

    void operator =( const Self& other )
    {
      this->m_Face = other.m_Face;
      this->m_Value = other.m_Value;
      this->m_Cluster = other.m_Cluster;
    }
    };

  typedef MinPriorityQueueElementWrapper<Item, MetricValueType,
                                         CellIdentifier>
    PriorityQueueItemType;
  typedef PriorityQueueContainer<PriorityQueueItemType,
                                 PriorityQueueItemType,
                                 MetricValueType, long>   PriorityQueueType;
  typedef typename PriorityQueueType::Pointer PriorityQueuePointer;
  typedef std::map<CellIdentifier, std::list<PriorityQueueItemType> >
    QueueMapType;
  typedef typename QueueMapType::iterator QueueMapIterator;

  itkSetConstObjectMacro( Mesh, MeshType );
  itkGetConstObjectMacro( Mesh, MeshType );

  void Reset()
  {
    m_FastMarchingComputed = false;
    m_ElementMap.clear();
    m_ClusterVector.clear();
    m_FaceClusterMap.clear();
    m_QueueMap.clear();
    m_ClusterBorderVector.clear();
    m_PriorityQueue->Clear();
    m_FaceProcessed.clear();
  }

  void SetSeedFaces( const SeedVectorType& iSeeds )
  {
    m_SeedFaces = iSeeds;
  }

  SeedVectorType GetSeedFaces() const
  {
    return m_SeedFaces;
  }

  CellIdentifier GetSeedFace( const size_t& iId ) const
  {
    assert( ( iId >= 0 ) && ( iId < m_SeedFaces.size() ) );
    return m_SeedFaces[iId];
  }

  PointType GetSeedCaracteristicPoint( const size_t& iId ) const
  {
    CellAutoPointer cell1;

    m_Mesh->GetCell( GetSeedFace( iId ), cell1 );
    PointIdIterator it1 = cell1->PointIdsBegin();

    PointType pt1[3];
    for( int k = 0; k < 3; it1++, k++ )
      {
      pt1[k] = m_Mesh->GetPoint( *it1 );
      }

    return TriangleType::ComputeGravityCenter( pt1[0], pt1[1], pt1[2] );
  }

  ClusterVectorType GetCluster()
  {
    assert( m_Mesh.IsNotNull() );
    assert( m_SeedFaces.size() > 1 );

    m_Metric.SetMesh( m_Mesh );

    if( !m_FastMarchingComputed )
      {
      ComputeFastMarching();
      }

    return m_ClusterVector;
  }

  ClusterType Evaluate( MeshType* iMesh,
                        const SeedVectorType& iSeeds,
                        const size_t& iId )
  {
    assert( iMesh != 0 );
    assert( iSeeds.size() > 1 );
    assert( ( iId >= 0 ) && ( iId < iSeeds.size() ) );

    m_Mesh = iMesh;
    m_Metric.SetMesh( m_Mesh );
    m_SeedFaces = iSeeds;

    if( !m_FastMarchingComputed )
      {
      ComputeFastMarching();
      }

    return m_ClusterVector[iId];
  }

  ClusterType Evaluate( const size_t& iId )
  {
    assert( m_Mesh.IsNotNull() );
    assert( m_SeedFaces.size() > 1 );
    assert( iId < m_SeedFaces.size() );

    m_Metric.SetMesh( m_Mesh );

    if( !m_FastMarchingComputed )
      {
      ComputeFastMarching();
      }

    return m_ClusterVector[iId];
  }

  ClusterType Evaluate( const SeedVectorType& iSeeds,
                        const size_t& iId )
  {
    assert( m_Mesh.IsNotNull() );
    assert( iSeeds.size() > 1 );
    assert( ( iId >= 0 ) && ( iId < iSeeds.size() ) );

    m_SeedFaces = iSeeds;
    m_Metric.SetMesh( m_Mesh );

    if( !m_FastMarchingComputed )
      {
      ComputeFastMarching();
      }

    return m_ClusterVector[iId];
  }

  std::vector<CellIdentifierListType> m_ClusterBorderVector;
protected:
  QuadEdgeMeshDualFastMarching() : m_Mesh( nullptr ),
    m_NumberOfUnprocessedElements( 1 ),
    m_FastMarchingComputed( false )
  {
    m_PriorityQueue = PriorityQueueType::New();
  }

  ~QuadEdgeMeshDualFastMarching()
  {
  }

  MeshConstPointer                 m_Mesh;
  MetricType                       m_Metric;
  SeedVectorType                   m_SeedFaces;
  std::map<CellIdentifier, bool>   m_FaceProcessed;
  QueueMapType                     m_ElementMap;
  ClusterVectorType                m_ClusterVector;
  std::map<CellIdentifier, size_t> m_FaceClusterMap;
//         std::vector< std::list< CellIdentifier> > m_ClusterBorderVector;
  PriorityQueuePointer m_PriorityQueue;
  CellIdentifier       m_NumberOfUnprocessedElements;
  bool                 m_FastMarchingComputed;
  QueueMapType         m_QueueMap;

  void InitFMM()
  {
    m_NumberOfUnprocessedElements = m_Mesh->GetNumberOfCells();

    CellsContainerConstPointer meshcells = m_Mesh->GetCells();

    CellsContainerConstIterator cell_it = meshcells->Begin();

    while( cell_it != meshcells->End() )
      {
      if( cell_it.Value()->GetNumberOfPoints() > 2 )
        {
        m_FaceProcessed[cell_it->Index()] = false;
        }
      else
        {
        m_FaceProcessed[cell_it->Index()] = true;
        m_NumberOfUnprocessedElements--;
        }
      ++cell_it;
      }

    m_ClusterVector.resize( m_SeedFaces.size() );
    m_ClusterBorderVector.resize( m_SeedFaces.size() );

    CellIdentifier id_cluster( 0 ), id_face( 0 );

    SeedVectorIterator it = m_SeedFaces.begin();

    while( it != m_SeedFaces.end() )
      {
      id_face = *it;

      m_FaceProcessed[id_face] = true;

      Item                  item( id_face, 0., id_cluster );
      PriorityQueueItemType qi( item, 0. );
      m_QueueMap[id_face].push_back( qi );

      // push the first cell (seed) into cluster
      m_ClusterVector[id_cluster].push_back( id_face );
      m_FaceClusterMap[id_face] = id_cluster;

      PushFMM( id_face, id_cluster, 0. );

      ++it;
      ++id_cluster;
      --m_NumberOfUnprocessedElements;
      }
  }

  void PushFMM( const CellIdentifier& iId_face,
                const CellIdentifier& iId_cluster,
                const MetricValueType& iValue )
  {
    PolygonType* poly =
      dynamic_cast<PolygonType *>( m_Mesh->GetCells(
                                     )->GetElement( iId_face ) );

    if( poly != nullptr )
      {
      QEType* edge = poly->GetEdgeRingEntry();

      if( edge != nullptr )
        {
        if( edge->GetLeft() != iId_face )
          {
          edge = edge->GetSym();
          }

        // ***
        QEType *        temp( edge );
        MetricValueType value( iValue );
        CellIdentifier  id_face2( 0 ), id_cluster2( 0 );

        do
          {
          id_face2 = temp->GetRight();

          if( id_face2 != MeshType::m_NoFace )
            {
            if( !m_FaceProcessed[id_face2] )
              {
              value += m_Metric( iId_face, id_face2 );

              Item                  item( id_face2, value, iId_cluster );
              PriorityQueueItemType qj( item, value );

              m_QueueMap[iId_face].push_back( qj );
              m_PriorityQueue->Push( qj );
              }
            else
              {
              id_cluster2 = m_FaceClusterMap[id_face2];
              if( id_cluster2 != iId_cluster )
                {
                m_ClusterBorderVector[iId_cluster].push_back( iId_face );
                m_ClusterBorderVector[id_cluster2].push_back( id_face2 );
                }
              }
            }
          temp = temp->GetLnext();
          }
        while( temp != edge );

//               }
        }
      }
  }

  void UpdateFMM()
  {
    MetricValueType       value( 0. );
    CellIdentifier        id_cluster( 0 ), id_face( 0 );
    PriorityQueueItemType qi;
    Item                  item;

    while( !m_PriorityQueue->Empty() )
      {
      qi = m_PriorityQueue->Peek();
      item = qi.m_Element;
      id_face = item.m_Face;
      id_cluster = item.m_Cluster;
      value = item.m_Value;

      if( !m_FaceProcessed[id_face] )
        {
        QueueMapIterator queue_it = m_QueueMap.find( id_face );

        if( queue_it != m_QueueMap.end() )
          {
          while( !queue_it->second.empty() )
            {
            m_PriorityQueue->DeleteElement( queue_it->second.front() );
            queue_it->second.pop_front();
            }

          m_QueueMap.erase( queue_it );
          }

        m_FaceProcessed[id_face] = true;
        // push more cell into the cluster
        m_ClusterVector[id_cluster].push_back( id_face );
        m_FaceClusterMap[id_face] = id_cluster;

        PushFMM( id_face, id_cluster, value );
        }
      m_PriorityQueue->Pop();
      --m_NumberOfUnprocessedElements;
      }
  }

  void ComputeFastMarching()
  {
    InitFMM();
    UpdateFMM();
    m_FastMarchingComputed = true;
  }

private:
  QuadEdgeMeshDualFastMarching( const Self & );
  void operator =( const Self & );
};
}

#endif
