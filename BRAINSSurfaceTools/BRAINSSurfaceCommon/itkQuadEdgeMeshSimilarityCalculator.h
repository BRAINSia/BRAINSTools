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
template <class TInputMesh1, class TInputMesh2 = TInputMesh1>
class QuadEdgeMeshSimilarityCalculator :
  public QuadEdgeMeshToQuadEdgeMeshFilter<TInputMesh1, TInputMesh2>
{
public:
  typedef QuadEdgeMeshSimilarityCalculator                           Self;
  typedef QuadEdgeMeshToQuadEdgeMeshFilter<TInputMesh1, TInputMesh2> Superclass;
  typedef SmartPointer<Self>                                         Pointer;
  typedef SmartPointer<const Self>                                   ConstPointer;

  /** Run-time type information (and related methods).   */
  itkTypeMacro( QuadEdgeMeshSimilarityCalculator, QuadEdgeMeshToQuadEdgeMeshFilter );

  /** New macro for creation of through a Smart Pointer   */
  itkNewMacro( Self );

  typedef TInputMesh1                                     InputMeshType1;
  typedef typename InputMeshType1::PointType              InputPointType1;
  typedef typename InputMeshType1::PixelType              InputPixelType1;
  typedef typename InputMeshType1::PointIdentifier        PointIdentifier1;
  typedef typename InputMeshType1::PointsContainer        InputPointsContainer1;
  typedef typename InputMeshType1::PointsContainerPointer InputPointsContainerPointer1;
  typedef typename InputMeshType1::PointDataContainer     InputPointDataContainer1;

  typedef TInputMesh2                                     InputMeshType2;
  typedef typename InputMeshType2::PointType              InputPointType2;
  typedef typename InputMeshType2::PixelType              InputPixelType2;
  typedef typename InputMeshType2::PointIdentifier        PointIdentifier2;
  typedef typename InputMeshType2::PointsContainer        InputPointsContainer2;
  typedef typename InputMeshType2::PointsContainerPointer InputPointsContainerPointer2;
  typedef typename InputMeshType2::PointDataContainer     InputPointDataContainer2;

  typedef typename InputMeshType1::CellType       CellType1;
  typedef typename InputMeshType1::CellTraits     CellTraits1;
  typedef typename InputMeshType1::CellsContainer CellsContainer1;
  typedef typename CellsContainer1::Iterator      CellsContainerIterator1;
  typedef typename CellsContainer1::ConstIterator CellsContainerConstIterator1;
  typedef typename CellTraits1::PointIdIterator   PointIdIterator1;

  typedef typename InputMeshType2::CellType       CellType2;
  typedef typename InputMeshType2::CellTraits     CellTraits2;
  typedef typename InputMeshType2::CellsContainer CellsContainer2;
  typedef typename CellsContainer2::Iterator      CellsContainerIterator2;
  typedef typename CellsContainer2::ConstIterator CellsContainerConstIterator2;
  typedef typename CellTraits2::PointIdIterator   PointIdIterator2;

  typedef TriangleHelper<InputPointType1>     TriangleType;
  typedef typename TriangleType::CoordRepType AreaType;

  /** Set/Get the mesh1 that has the labels. */
  void SetInputMesh1( const InputMeshType1 * mesh1 );

  const InputMeshType1 * GetInputMesh1( void ) const;

  /** Set/Get the mesh2 that has the labels. */
  void SetInputMesh2( const InputMeshType2 * mesh2 );

  const InputMeshType2 * GetInputMesh2( void ) const;

  /** Set the label value we want to calculate. */
  itkSetMacro( LabelValue, InputPixelType1 );
  /** Get the label value it is calculating. */
  itkGetMacro( LabelValue, InputPixelType1 );

  /** Get the calculated Dice Value. */
  itkGetMacro( Dice, double );
  /** Get the calculated Jaccard Value. */
  itkGetMacro( Jaccard, double );

  void Compute();

protected:
  QuadEdgeMeshSimilarityCalculator();
  ~QuadEdgeMeshSimilarityCalculator();
private:

  ITK_DISALLOW_COPY_AND_ASSIGN(QuadEdgeMeshSimilarityCalculator);

  InputPixelType1 m_LabelValue;

  double m_Dice;
  double m_Jaccard;
};
}

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkQuadEdgeMeshSimilarityCalculator.hxx"
#endif

#endif
