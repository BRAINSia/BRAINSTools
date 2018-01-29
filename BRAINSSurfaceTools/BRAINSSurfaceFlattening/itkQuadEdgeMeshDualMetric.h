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
#ifndef __itkQuadEdgeMeshDualMetric_h
#define __itkQuadEdgeMeshDualMetric_h

#include <itkTriangleHelper.h>

namespace itk
{
/**
*  \class QuadEdgeMeshDualMetric
*  \brief Functor which computes a distance between two dual cells.
* */
template <class TMesh>
class QuadEdgeMeshDualMetric
{
public:
  typedef TMesh                              MeshType;
  typedef typename MeshType::ConstPointer    MeshConstPointer;
  typedef typename MeshType::CellIdentifier  CellIdentifier;
  typedef typename MeshType::CellAutoPointer CellAutoPointer;
  typedef typename MeshType::PointIdIterator PointIdIterator;
  typedef typename MeshType::PointType       PointType;
  typedef typename PointType::RealType       ValueType;

  QuadEdgeMeshDualMetric() : m_Mesh( nullptr )
  {
  }

  QuadEdgeMeshDualMetric( MeshType* iMesh ) : m_Mesh( iMesh )
  {
  }

  virtual ~QuadEdgeMeshDualMetric()
  {
  }

  void SetMesh( const MeshType* iMesh )
  {
    m_Mesh = iMesh;
  }

  virtual ValueType operator ()( constexpr CellIdentifier& iFace1, constexpr CellIdentifier& iFace2 ) const = 0;

protected:
  MeshConstPointer m_Mesh;
};

/** \class QuadEdgeMeshDualSquaredEuclideanMetric
* \brief Functor which computes the distance between 2 faces as the squared
* Euclidean distance between their center of mass.
* */
template <class TMesh>
class QuadEdgeMeshDualSquaredEuclideanMetric : public
  QuadEdgeMeshDualMetric<TMesh>
{
public:
  typedef QuadEdgeMeshDualMetric<TMesh>          Superclass;
  typedef QuadEdgeMeshDualSquaredEuclideanMetric Self;

  typedef typename Superclass::MeshType         MeshType;
  typedef typename Superclass::MeshConstPointer MeshConstPointer;
  typedef typename Superclass::CellIdentifier   CellIdentifier;
  typedef typename Superclass::CellAutoPointer  CellAutoPointer;
  typedef typename Superclass::PointIdIterator  PointIdIterator;
  typedef typename Superclass::PointType        PointType;
  typedef typename Superclass::ValueType        ValueType;

  typedef TriangleHelper<PointType> TriangleType;

  QuadEdgeMeshDualSquaredEuclideanMetric() : Superclass()
  {
  }

  QuadEdgeMeshDualSquaredEuclideanMetric( MeshType* iMesh ) :
    Superclass( iMesh )
  {
  }

  ~QuadEdgeMeshDualSquaredEuclideanMetric()
  {
  }

  ValueType operator ()( const CellIdentifier& iFace1,
                         const CellIdentifier& iFace2 ) const
  {
    CellAutoPointer cell1, cell2;

    this->m_Mesh->GetCell( iFace1, cell1 );
    this->m_Mesh->GetCell( iFace2, cell2 );

    PointIdIterator it1 = cell1->PointIdsBegin();
    PointIdIterator it2 = cell2->PointIdsBegin();

    PointType pt1[3], pt2[3];
    for( int k = 0; k < 3; ++it1, ++it2, ++k )
      {
      pt1[k] = this->m_Mesh->GetPoint( *it1 );
      pt2[k] = this->m_Mesh->GetPoint( *it2 );
      }
    PointType center1 =
      TriangleType::ComputeGravityCenter( pt1[0], pt1[1], pt1[2] );

    PointType center2 =
      TriangleType::ComputeGravityCenter( pt2[0], pt2[1], pt2[2] );

    return center1.SquaredEuclideanDistanceTo( center2 );
  }
};

/** \class QuadEdgeMeshDualSquaredEuclideanWithAreaWeightMetric
* \brief Functors which computes the distance between two faces as the
* squared Euclidean distance between their center of mass, weighted by the
* area of the first face.
* */
template <class TMesh>
class QuadEdgeMeshDualSquaredEuclideanWithAreaWeightMetric
  : public QuadEdgeMeshDualMetric<TMesh>
{
public:
  typedef QuadEdgeMeshDualMetric<TMesh>                        Superclass;
  typedef QuadEdgeMeshDualSquaredEuclideanWithAreaWeightMetric Self;

  typedef typename Superclass::MeshType         MeshType;
  typedef typename Superclass::MeshConstPointer MeshConstPointer;
  typedef typename Superclass::CellIdentifier   CellIdentifier;
  typedef typename Superclass::CellAutoPointer  CellAutoPointer;
  typedef typename Superclass::PointIdIterator  PointIdIterator;
  typedef typename Superclass::PointType        PointType;
  typedef typename Superclass::ValueType        ValueType;

  typedef TriangleHelper<PointType> TriangleType;

  QuadEdgeMeshDualSquaredEuclideanWithAreaWeightMetric() : Superclass()
  {
  }

  QuadEdgeMeshDualSquaredEuclideanWithAreaWeightMetric( MeshType* iMesh ) :
    Superclass( iMesh )
  {
  }

  ~QuadEdgeMeshDualSquaredEuclideanWithAreaWeightMetric()
  {
  }

  ValueType operator ()( const CellIdentifier& iOrg,
                         const CellIdentifier& iDest ) const
  {
    CellAutoPointer cell1, cell2;

    this->m_Mesh->GetCell( iOrg, cell1 );
    this->m_Mesh->GetCell( iDest, cell2 );

    PointIdIterator it1 = cell1->PointIdsBegin();
    PointIdIterator it2 = cell2->PointIdsBegin();

    PointType pt1[3], pt2[3];
    for( int k = 0; k < 3; ++it1, ++it2, ++k )
      {
      pt1[k] = this->m_Mesh->GetPoint( *it1 );
      pt2[k] = this->m_Mesh->GetPoint( *it2 );
      }
    typename TriangleType::Pointer t1 = TriangleType::New();
    t1->SetPoints( pt1[0], pt1[1], pt1[2] );

    typename TriangleType::Pointer t2 = TriangleType::New();
    t2->SetPoints( pt2[0], pt2[1], pt2[2] );

    PointType center1 = t1->ComputeGravityCenter();
    PointType center2 = t2->ComputeGravityCenter();

    ValueType d_2 = center1.SquaredEuclideanDistanceTo( center2 );
    ValueType Area = t1->ComputeArea();

    return Area * d_2;
  }
};
}

#endif
