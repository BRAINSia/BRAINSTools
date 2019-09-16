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
#ifndef __itkQuadEdgeMeshBoundarySmoothFilter_h
#define __itkQuadEdgeMeshBoundarySmoothFilter_h

#include <itkQuadEdgeMeshToQuadEdgeMeshFilter.h>

namespace itk
{
/**
 * \class QuadEdgeMeshBoundarySmoothFilter
 * \brief This filter takes as inout two surfaces that represent
 * a single genus 0 surface. The smoothing will generate a smooth
 * boundary between the two surfaces. The filter produces two outputs
 * that are the resulting smoothed surfaces.
 *
 * The filter is supposed to be used after itkQuadEdgeMeshSplitFilter
 * and before mapping the surface to sphere.

 * itkQuadEdgeMeshSplitFilter takes one surface as input and
 * split it into two surfaces. There are "teeth" left on boundaries
 * of splitted surfaces because of the triangulars on the surface.
 *
 * The smoothing method is repeated for the user specified number
 * of iterations or until there is no triangles to move between
 * surfaces.
 *
 * Note: the total number of triangles on both inputs have to be
 * the same after each iteration. Missing cell or extra cell to
 * the whole surface is not allowed.
 *
 * The input mesh must share the same mesh type, but the output
 * may differ.
 * \ingroup MeshFilters
 *
 * \sa itkQuadEdgeMeshSplitFilter
 *
 */
template <typename TInputMesh, typename TOutputMesh>
class QuadEdgeMeshBoundarySmoothFilter : public QuadEdgeMeshToQuadEdgeMeshFilter<TInputMesh, TOutputMesh>
{
public:
  using Self = QuadEdgeMeshBoundarySmoothFilter;
  using Pointer = SmartPointer<Self>;
  using ConstPointer = SmartPointer<const Self>;
  using Superclass = QuadEdgeMeshToQuadEdgeMeshFilter<TInputMesh, TOutputMesh>;

  /** Input types. */
  using InputMeshType = TInputMesh;
  using InputMeshPointer = typename InputMeshType::Pointer;
  using InputMeshConstPointer = typename InputMeshType::ConstPointer;
  using InputCoordRepType = typename InputMeshType::CoordRepType;
  using InputPointType = typename InputMeshType::PointType;
  using InputPointIdentifier = typename InputMeshType::PointIdentifier;
  using InputQEType = typename InputMeshType::QEType;
  using InputQEPrimal = typename InputMeshType::QEPrimal;
  using InputVectorType = typename InputMeshType::VectorType;
  using InputCellType = typename InputMeshType::CellType;
  using InputCellIdentifier = typename InputMeshType::CellIdentifier;

  typedef typename InputMeshType::PointsContainerConstPointer  InputPointsContainerConstPointer;
  typedef typename InputMeshType::PointsContainerConstIterator InputPointsContainerConstIterator;
  typedef typename InputMeshType::CellsContainerConstPointer   InputCellsContainerConstPointer;
  typedef typename InputMeshType::CellsContainerConstIterator  InputCellsContainerConstIterator;

  using InputEdgeCellType = typename InputMeshType::EdgeCellType;
  using InputPolygonCellType = typename InputMeshType::PolygonCellType;
  using InputPointIdList = typename InputMeshType::PointIdList;
  using InputCellTraits = typename InputMeshType::CellTraits;
  typedef typename InputCellTraits::PointIdInternalIterator InputPointsIdInternalIterator;
  //    using InputQEIterator = typename InputQEPrimal::IteratorGeom;

  /** Output types. */
  using OutputMeshType = TOutputMesh;
  using OutputMeshPointer = typename OutputMeshType::Pointer;
  using OutputMeshConstPointer = typename OutputMeshType::ConstPointer;
  using OutputCoordRepType = typename OutputMeshType::CoordRepType;
  using OutputPointType = typename OutputMeshType::PointType;
  using OutputPointIdentifier = typename OutputMeshType::PointIdentifier;
  using OutputPointIdListType = typename OutputMeshType::PointIdList;

  using EdgeListPointerType = typename OutputMeshType::EdgeListPointerType;
  using OutputQEPrimal = typename OutputMeshType::QEPrimal;
  using OutputQEIterator = typename OutputQEPrimal::IteratorGeom;
  using OutputPolygonCellType = typename OutputMeshType::PolygonCellType;
  typedef typename OutputMeshType::PointsContainerPointer      OutputPointsContainerPointer;
  typedef typename OutputMeshType::PointsContainerIterator     OutputPointsContainerIterator;
  typedef typename OutputMeshType::CellsContainerConstPointer  OutputCellsContainerConstPointer;
  typedef typename OutputMeshType::CellsContainerConstIterator OutputCellsContainerConstIterator;

public:
  itkNewMacro(Self);
  itkTypeMacro(QuadEdgeMeshBoundarySmoothFilter, QuadEdgeMeshToQuadEdgeMeshFilter);

  /** Set/Get the first input mesh */
  void
  SetInputMesh1(const InputMeshType * mesh1);

  const InputMeshType *
  GetInputMesh1(void) const;

  /** Set/Get the second input mesh */
  void
  SetInputMesh2(const InputMeshType * mesh2);

  const InputMeshType *
  GetInputMesh2(void) const;

  /** Get the first smoothed hemisphere */
  OutputMeshType *
  GetOutputMesh1(void);

  /** Get the second smoothed hemisphere */
  OutputMeshType *
  GetOutputMesh2(void);

  /** Set/Get the number of iterations. */
  itkSetMacro(Iterations, int);
  itkGetMacro(Iterations, int);

protected:
  QuadEdgeMeshBoundarySmoothFilter();
  ~QuadEdgeMeshBoundarySmoothFilter();

  void
  CopyInputMeshesToOutputMeshes();

  int
  AdjustBoundary(OutputMeshType * deleteMesh, OutputMeshType * addMesh);

  virtual void
  GenerateData() override;

private:
  QuadEdgeMeshBoundarySmoothFilter(const Self &);
  void
  operator=(const Self &);

  int m_Iterations;
};
} // namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#  include "itkQuadEdgeMeshBoundarySmoothFilter.hxx"
#endif

#endif
