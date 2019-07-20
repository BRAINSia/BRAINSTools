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
#ifndef __itkQuadEdgeMeshSimilarityCalculator_h
#define __itkQuadEdgeMeshSimilarityCalculator_h

#include "itkQuadEdgeMeshToQuadEdgeMeshFilter.h"

#include "itkTriangleHelper.h"

namespace itk
{
/**
 * \class QuadEdgeMeshSimilarityCalculator
 * \brief The calculator takes two input meshes and give
 *  Dice or Jaccard Similarity by users' demand.
 *
 * QuadEdgeMeshSimilarityCalculator requires two input meshes
 * have the same topologies (points, cells and edges are all the same).
 *
 * Labels are given as integer scalars on meshes.
 *
 * Similarities are calculated according to the area of labels.
 *
 * The area of a point's label in a triangle is calculated as 1/3 of the
 * whole triangle.
 *
 * Area of a particular label is the sum of area of that label in all triangles.
 *
 * \ingroup MeshFilters
 *
 */
template < typename TInputMesh1, typename TInputMesh2 = TInputMesh1 >
class QuadEdgeMeshSimilarityCalculator : public QuadEdgeMeshToQuadEdgeMeshFilter< TInputMesh1, TInputMesh2 >
{
public:
  ITK_DISALLOW_COPY_AND_ASSIGN( QuadEdgeMeshSimilarityCalculator );

  using Self = QuadEdgeMeshSimilarityCalculator;
  using Superclass = QuadEdgeMeshToQuadEdgeMeshFilter< TInputMesh1, TInputMesh2 >;
  using Pointer = SmartPointer< Self >;
  using ConstPointer = SmartPointer< const Self >;

  /** Run-time type information (and related methods).   */
  itkTypeMacro( QuadEdgeMeshSimilarityCalculator, QuadEdgeMeshToQuadEdgeMeshFilter );

  /** New macro for creation of through a Smart Pointer   */
  itkNewMacro( Self );

  using InputMeshType1 = TInputMesh1;
  using InputPointType1 = typename InputMeshType1::PointType;
  using InputPixelType1 = typename InputMeshType1::PixelType;
  using PointIdentifier1 = typename InputMeshType1::PointIdentifier;
  using InputPointsContainer1 = typename InputMeshType1::PointsContainer;
  using InputPointsContainerPointer1 = typename InputMeshType1::PointsContainerPointer;
  using InputPointDataContainer1 = typename InputMeshType1::PointDataContainer;

  using InputMeshType2 = TInputMesh2;
  using InputPointType2 = typename InputMeshType2::PointType;
  using InputPixelType2 = typename InputMeshType2::PixelType;
  using PointIdentifier2 = typename InputMeshType2::PointIdentifier;
  using InputPointsContainer2 = typename InputMeshType2::PointsContainer;
  using InputPointsContainerPointer2 = typename InputMeshType2::PointsContainerPointer;
  using InputPointDataContainer2 = typename InputMeshType2::PointDataContainer;

  using CellType1 = typename InputMeshType1::CellType;
  using CellTraits1 = typename InputMeshType1::CellTraits;
  using CellsContainer1 = typename InputMeshType1::CellsContainer;
  using CellsContainerIterator1 = typename CellsContainer1::Iterator;
  using CellsContainerConstIterator1 = typename CellsContainer1::ConstIterator;
  using PointIdIterator1 = typename CellTraits1::PointIdIterator;

  using CellType2 = typename InputMeshType2::CellType;
  using CellTraits2 = typename InputMeshType2::CellTraits;
  using CellsContainer2 = typename InputMeshType2::CellsContainer;
  using CellsContainerIterator2 = typename CellsContainer2::Iterator;
  using CellsContainerConstIterator2 = typename CellsContainer2::ConstIterator;
  using PointIdIterator2 = typename CellTraits2::PointIdIterator;

  using TriangleType = TriangleHelper< InputPointType1 >;
  using AreaType = typename TriangleType::CoordRepType;

  /** Set/Get the mesh1 that has the labels. */
  void
  SetInputMesh1( const InputMeshType1 * mesh1 );

  const InputMeshType1 *
  GetInputMesh1( void ) const;

  /** Set/Get the mesh2 that has the labels. */
  void
  SetInputMesh2( const InputMeshType2 * mesh2 );

  const InputMeshType2 *
  GetInputMesh2( void ) const;

  /** Set the label value we want to calculate. */
  itkSetMacro( LabelValue, InputPixelType1 );
  /** Get the label value it is calculating. */
  itkGetMacro( LabelValue, InputPixelType1 );

  /** Get the calculated Dice Value. */
  itkGetMacro( Dice, double );
  /** Get the calculated Jaccard Value. */
  itkGetMacro( Jaccard, double );

  void
  Compute();

protected:
  QuadEdgeMeshSimilarityCalculator();
  ~QuadEdgeMeshSimilarityCalculator();

private:
  InputPixelType1 m_LabelValue;

  double m_Dice;
  double m_Jaccard;
};
} // namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#  include "itkQuadEdgeMeshSimilarityCalculator.hxx"
#endif

#endif
