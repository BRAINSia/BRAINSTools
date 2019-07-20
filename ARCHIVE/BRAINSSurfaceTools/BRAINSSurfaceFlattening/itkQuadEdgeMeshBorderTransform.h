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
#ifndef __itkQuadEdgeMeshBorderTransform_h
#define __itkQuadEdgeMeshBorderTransform_h

#include <itkQuadEdgeMesh.h>
#include <itkQuadEdgeMeshToQuadEdgeMeshFilter.h>
#include <itkQuadEdgeMeshBoundaryEdgesMeshFunction.h>

namespace itk
{
/**
 * \class QuadEdgeMeshBorderTransform
 * \brief Transform the mandatoryly unique border of an \ref itkQE::Mesh
 * into either a circle (conformal) or a square (arclenght-wise).
 *
 * To Write.
 */
template < typename TInputMesh, typename TOutputMesh >
class QuadEdgeMeshBorderTransform : public QuadEdgeMeshToQuadEdgeMeshFilter< TInputMesh, TOutputMesh >
{
public:
  /** Basic types. */
  using Self = QuadEdgeMeshBorderTransform;
  using Superclass = QuadEdgeMeshToQuadEdgeMeshFilter< TInputMesh, TOutputMesh >;
  using Pointer = SmartPointer< Self >;
  using ConstPointer = SmartPointer< const Self >;

  using InputMeshType = TInputMesh;
  using InputMeshConstPointer = typename InputMeshType::ConstPointer;
  using InputCoordRepType = typename InputMeshType::CoordRepType;
  using InputPointType = typename InputMeshType::PointType;
  using InputTraits = typename InputMeshType::Traits;
  using InputPointIdentifier = typename InputMeshType::PointIdentifier;
  using InputQEType = typename InputMeshType::QEType;
  using InputIteratorGeom = typename InputQEType::IteratorGeom;
  using InputVectorType = typename InputMeshType::VectorType;
  using InputEdgeListType = typename InputMeshType::EdgeListType;
  using InputEdgeListPointerType = typename InputMeshType::EdgeListPointerType;
  using InputEdgeListIterator = typename InputEdgeListType::iterator;
  using InputEdgeCellType = typename InputMeshType::EdgeCellType;
  using InputPolygonCellType = typename InputMeshType::PolygonCellType;
  using InputPointIdList = typename InputMeshType::PointIdList;
  using InputPointsContainer = typename InputMeshType::PointsContainer;
  typedef typename InputMeshType::PointsContainerConstIterator InputPointsContainerConstIterator;
  typedef typename InputMeshType::CellsContainerConstIterator  InputCellsContainerConstIterator;

  using OutputMeshType = TOutputMesh;
  using OutputMeshPointer = typename OutputMeshType::Pointer;
  using OutputCoordRepType = typename OutputMeshType::CoordRepType;
  using OutputPointType = typename OutputMeshType::PointType;
  using OutputTraits = typename OutputMeshType::Traits;
  using OutputPointIdentifier = typename OutputMeshType::PointIdentifier;
  using OutputQEType = typename OutputMeshType::QEType;
  using OutputVectorType = typename OutputMeshType::VectorType;
  using OutputEdgeListType = typename OutputMeshType::EdgeListType;
  using OutputEdgeCellType = typename OutputMeshType::EdgeCellType;
  using OutputPolygonCellType = typename OutputMeshType::PolygonCellType;
  using OutputPointIdList = typename OutputMeshType::PointIdList;
  using OutputPointsContainer = typename OutputMeshType::PointsContainer;
  typedef typename OutputMeshType::PointsContainerConstIterator OutputPointsContainerConstIterator;
  typedef typename OutputMeshType::CellsContainerConstIterator  OutputCellsContainerConstIterator;

  itkNewMacro( Self );
  itkTypeMacro( QuadEdgeMeshBorderTransform, QuadEdgeMeshToQuadEdgeMeshFilter );
  static constexpr unsigned int PointDimension = InputTraits::PointDimension;

  using InputVectorPointType = std::vector< InputPointType >;
  using MapPointIdentifier = std::map< InputPointIdentifier, OutputPointIdentifier >;
  using MapPointIdentifierIterator = typename MapPointIdentifier::iterator;

  using BoundaryRepresentativeEdgesType = QuadEdgeMeshBoundaryEdgesMeshFunction< InputMeshType >;
  using BoundaryRepresentativeEdgesPointer = typename BoundaryRepresentativeEdgesType::Pointer;

public:
  enum BorderTransformType
  {
    SQUARE_BORDER_TRANSFORM = 0,
    DISK_BORDER_TRANSFORM
  };

  itkSetMacro( TransformType, BorderTransformType );
  itkGetConstMacro( TransformType, BorderTransformType );

  itkSetMacro( Radius, InputCoordRepType );
  itkGetConstMacro( Radius, InputCoordRepType );

  void
  ComputeTransform();

  MapPointIdentifier
  GetBoundaryPtMap();

  InputVectorPointType
  GetBorder();

protected:
  QuadEdgeMeshBorderTransform();
  ~QuadEdgeMeshBorderTransform(){};

  BorderTransformType m_TransformType;

  InputCoordRepType    m_Radius;
  InputVectorPointType m_Border;

  MapPointIdentifier m_BoundaryPtMap;

  void
  GenerateData() override;

  void
  ComputeBoundary();

  InputEdgeListIterator
  ComputeLongestBorder();

  InputEdgeListIterator
  ComputeLargestBorder();

  void
  DiskTransform();

  InputPointType
  GetMeshBarycentre();

  InputCoordRepType
  RadiusMaxSquare();

  void
  ArcLengthSquareTransform();

private:
  /** Not implemented */
  QuadEdgeMeshBorderTransform( const Self & );

  /** Not implemented */
  void
  operator=( const Self & );
};
} // namespace itk
#include "itkQuadEdgeMeshBorderTransform.hxx"

#endif
