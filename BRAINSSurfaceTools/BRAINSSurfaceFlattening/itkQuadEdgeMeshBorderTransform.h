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
template <class TInputMesh, class TOutputMesh>
class QuadEdgeMeshBorderTransform :
  public QuadEdgeMeshToQuadEdgeMeshFilter<TInputMesh, TOutputMesh>
{
public:
  /** Basic types. */
  typedef QuadEdgeMeshBorderTransform Self;
  typedef QuadEdgeMeshToQuadEdgeMeshFilter<TInputMesh,
                                           TOutputMesh>                                       Superclass;
  typedef SmartPointer<Self>       Pointer;
  typedef SmartPointer<const Self> ConstPointer;

  typedef TInputMesh                                  InputMeshType;
  typedef typename InputMeshType::ConstPointer        InputMeshConstPointer;
  typedef typename InputMeshType::CoordRepType        InputCoordRepType;
  typedef typename InputMeshType::PointType           InputPointType;
  typedef typename InputMeshType::Traits              InputTraits;
  typedef typename InputMeshType::PointIdentifier     InputPointIdentifier;
  typedef typename InputMeshType::QEType              InputQEType;
  typedef typename InputQEType::IteratorGeom          InputIteratorGeom;
  typedef typename InputMeshType::VectorType          InputVectorType;
  typedef typename InputMeshType::EdgeListType        InputEdgeListType;
  typedef typename InputMeshType::EdgeListPointerType InputEdgeListPointerType;
  typedef typename InputEdgeListType::iterator        InputEdgeListIterator;
  typedef typename InputMeshType::EdgeCellType        InputEdgeCellType;
  typedef typename InputMeshType::PolygonCellType     InputPolygonCellType;
  typedef typename InputMeshType::PointIdList         InputPointIdList;
  typedef typename InputMeshType::PointsContainer     InputPointsContainer;
  typedef typename InputMeshType::PointsContainerConstIterator
    InputPointsContainerConstIterator;
  typedef typename InputMeshType::CellsContainerConstIterator
    InputCellsContainerConstIterator;

  typedef TOutputMesh                              OutputMeshType;
  typedef typename OutputMeshType::Pointer         OutputMeshPointer;
  typedef typename OutputMeshType::CoordRepType    OutputCoordRepType;
  typedef typename OutputMeshType::PointType       OutputPointType;
  typedef typename OutputMeshType::Traits          OutputTraits;
  typedef typename OutputMeshType::PointIdentifier OutputPointIdentifier;
  typedef typename OutputMeshType::QEType          OutputQEType;
  typedef typename OutputMeshType::VectorType      OutputVectorType;
  typedef typename OutputMeshType::EdgeListType    OutputEdgeListType;
  typedef typename OutputMeshType::EdgeCellType    OutputEdgeCellType;
  typedef typename OutputMeshType::PolygonCellType OutputPolygonCellType;
  typedef typename OutputMeshType::PointIdList     OutputPointIdList;
  typedef typename OutputMeshType::PointsContainer OutputPointsContainer;
  typedef typename OutputMeshType::PointsContainerConstIterator
    OutputPointsContainerConstIterator;
  typedef typename OutputMeshType::CellsContainerConstIterator
    OutputCellsContainerConstIterator;

  itkNewMacro( Self );
  itkTypeMacro( QuadEdgeMeshBorderTransform, QuadEdgeMeshToQuadEdgeMeshFilter );
  itkStaticConstMacro( PointDimension, unsigned int,
                       InputTraits::PointDimension );

  typedef std::vector<InputPointType>                           InputVectorPointType;
  typedef std::map<InputPointIdentifier, OutputPointIdentifier> MapPointIdentifier;
  typedef typename MapPointIdentifier::iterator                 MapPointIdentifierIterator;

  typedef QuadEdgeMeshBoundaryEdgesMeshFunction<InputMeshType> BoundaryRepresentativeEdgesType;
  typedef typename BoundaryRepresentativeEdgesType::Pointer    BoundaryRepresentativeEdgesPointer;
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

  void ComputeTransform();

  MapPointIdentifier GetBoundaryPtMap();

  InputVectorPointType GetBorder();

protected:
  QuadEdgeMeshBorderTransform();
  ~QuadEdgeMeshBorderTransform()
  {
  };

  BorderTransformType m_TransformType;

  InputCoordRepType    m_Radius;
  InputVectorPointType m_Border;

  MapPointIdentifier m_BoundaryPtMap;

  void GenerateData() override;

  void ComputeBoundary();

  InputEdgeListIterator ComputeLongestBorder();

  InputEdgeListIterator ComputeLargestBorder();

  void DiskTransform();

  InputPointType GetMeshBarycentre();

  InputCoordRepType RadiusMaxSquare();

  void ArcLengthSquareTransform();

private:
  /** Not implemented */
  QuadEdgeMeshBorderTransform( const Self & );

  /** Not implemented */
  void operator =( const Self & );
};
}
#include "itkQuadEdgeMeshBorderTransform.hxx"

#endif
