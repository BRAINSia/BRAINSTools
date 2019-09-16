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

#ifndef __itkQuadEdgeMeshSplitFilter_h
#define __itkQuadEdgeMeshSplitFilter_h

#include <itkQuadEdgeMeshToQuadEdgeMeshFilter.h>
#include <itkQuadEdgeMeshPolygonCell.h>
#include <itkTriangleHelper.h>

#include "itkQuadEdgeMeshDualFastMarching.h"

namespace itk
{
template <typename TInputMesh, typename TOutputMesh>
class QuadEdgeMeshSplitFilter : public QuadEdgeMeshToQuadEdgeMeshFilter<TInputMesh, TOutputMesh>
{
public:
  using Self = QuadEdgeMeshSplitFilter;
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
  using SeedVectorType = std::vector<InputCellIdentifier>;

  using InputPolygonType = QuadEdgeMeshPolygonCell<InputCellType>;
  using InputPolygonAutoPointer = typename InputPolygonType::SelfAutoPointer;

  using PolygonSetType = std::list<InputPolygonType *>;

  using FaceAreaMapType = std::map<InputPolygonType *, InputCoordRepType>;

  using TriangleType = TriangleHelper<InputPointType>;

  using FMMType = QuadEdgeMeshDualFastMarching<InputMeshType>;
  using FMMPointer = typename FMMType::Pointer;
  using FMMClusterType = typename FMMType::ClusterType;
  using FMMClusterIterator = typename FMMType::ClusterIterator;
  using FMMSeedVectorType = typename FMMType::SeedVectorType;

  /** Output types. */
  using OutputMeshType = TOutputMesh;
  using OutputMeshPointer = typename OutputMeshType::Pointer;
  using OutputMeshConstPointer = typename OutputMeshType::ConstPointer;
  using OutputCoordRepType = typename OutputMeshType::CoordRepType;
  using OutputPointType = typename OutputMeshType::PointType;
  using OutputPointIdentifier = typename OutputMeshType::PointIdentifier;
  using OutputPointIdList = typename OutputMeshType::PointIdList;
  using OutputQEType = typename OutputMeshType::QEType;
  using OutputVectorType = typename OutputMeshType::VectorType;
  using OutputCellType = typename OutputMeshType::CellType;
  //    using OutputQEIterator = typename OutputQEPrimal::IteratorGeom;
  typedef typename OutputMeshType::PointsContainerPointer  OutputPointsContainerPointer;
  typedef typename OutputMeshType::PointsContainerIterator OutputPointsContainerIterator;

  using OutputPolygonType = QuadEdgeMeshPolygonCell<OutputCellType>;
  using OutputPolygonAutoPointer = typename OutputPolygonType::SelfAutoPointer;

public:
  itkNewMacro(Self);
  itkTypeMacro(QuadEdgeMeshSplitFilter, QuadEdgeMeshToQuadEdgeMeshFilter);

  void
  SetSeedFaces(const SeedVectorType & iSeeds)
  {
    m_SeedFaces = iSeeds;
  }

protected:
  QuadEdgeMeshSplitFilter();
  ~QuadEdgeMeshSplitFilter();

  InputCellIdentifier m_StartCellId;
  InputCoordRepType   m_Area;
  FaceAreaMapType     m_FaceAreaMap;
  PolygonSetType      m_Faces;
  PolygonSetType      m_ProcessedFaces;
  SeedVectorType      m_SeedFaces;

  virtual void
  GenerateData() override;

  OutputPointIdList
  AddFacePointsToOutputMesh(OutputMeshType * iMesh, InputPolygonType * iPoly);

private:
  QuadEdgeMeshSplitFilter(const Self &);
  void
  operator=(const Self &);
};
} // namespace itk

#include "itkQuadEdgeMeshSplitFilter.hxx"
#endif
