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
template <class TInputMesh, class TOutputMesh>
class QuadEdgeMeshSplitFilter :
  public QuadEdgeMeshToQuadEdgeMeshFilter<TInputMesh, TOutputMesh>
{
public:
  typedef QuadEdgeMeshSplitFilter                                   Self;
  typedef SmartPointer<Self>                                        Pointer;
  typedef SmartPointer<const Self>                                  ConstPointer;
  typedef QuadEdgeMeshToQuadEdgeMeshFilter<TInputMesh, TOutputMesh> Superclass;

  /** Input types. */
  typedef TInputMesh                              InputMeshType;
  typedef typename InputMeshType::Pointer         InputMeshPointer;
  typedef typename InputMeshType::ConstPointer    InputMeshConstPointer;
  typedef typename InputMeshType::CoordRepType    InputCoordRepType;
  typedef typename InputMeshType::PointType       InputPointType;
  typedef typename InputMeshType::PointIdentifier InputPointIdentifier;
  typedef typename InputMeshType::QEType          InputQEType;
  typedef typename InputMeshType::VectorType      InputVectorType;
  typedef typename InputMeshType::CellType        InputCellType;
  typedef typename InputMeshType::CellIdentifier  InputCellIdentifier;

  typedef typename InputMeshType::PointsContainerConstPointer
    InputPointsContainerConstPointer;
  typedef typename InputMeshType::PointsContainerConstIterator
    InputPointsContainerConstIterator;
  typedef typename InputMeshType::CellsContainerConstPointer
    InputCellsContainerConstPointer;
  typedef typename InputMeshType::CellsContainerConstIterator
    InputCellsContainerConstIterator;

  typedef typename InputMeshType::EdgeCellType    InputEdgeCellType;
  typedef typename InputMeshType::PolygonCellType InputPolygonCellType;
  typedef typename InputMeshType::PointIdList     InputPointIdList;
  typedef typename InputMeshType::CellTraits      InputCellTraits;
  typedef typename InputCellTraits::PointIdInternalIterator
    InputPointsIdInternalIterator;
//    typedef typename InputQEPrimal::IteratorGeom    InputQEIterator;
  typedef std::vector<InputCellIdentifier> SeedVectorType;

  typedef QuadEdgeMeshPolygonCell<InputCellType>     InputPolygonType;
  typedef typename InputPolygonType::SelfAutoPointer InputPolygonAutoPointer;

  typedef std::list<InputPolygonType *> PolygonSetType;

  typedef std::map<InputPolygonType *, InputCoordRepType>
    FaceAreaMapType;

  typedef TriangleHelper<InputPointType> TriangleType;

  typedef QuadEdgeMeshDualFastMarching<InputMeshType> FMMType;
  typedef typename FMMType::Pointer                   FMMPointer;
  typedef typename FMMType::ClusterType               FMMClusterType;
  typedef typename FMMType::ClusterIterator           FMMClusterIterator;
  typedef typename FMMType::SeedVectorType            FMMSeedVectorType;

  /** Output types. */
  typedef TOutputMesh                              OutputMeshType;
  typedef typename OutputMeshType::Pointer         OutputMeshPointer;
  typedef typename OutputMeshType::ConstPointer    OutputMeshConstPointer;
  typedef typename OutputMeshType::CoordRepType    OutputCoordRepType;
  typedef typename OutputMeshType::PointType       OutputPointType;
  typedef typename OutputMeshType::PointIdentifier OutputPointIdentifier;
  typedef typename OutputMeshType::PointIdList     OutputPointIdList;
  typedef typename OutputMeshType::QEType          OutputQEType;
  typedef typename OutputMeshType::VectorType      OutputVectorType;
  typedef typename OutputMeshType::CellType        OutputCellType;
//    typedef typename OutputQEPrimal::IteratorGeom     OutputQEIterator;
  typedef typename OutputMeshType::PointsContainerPointer
    OutputPointsContainerPointer;
  typedef typename OutputMeshType::PointsContainerIterator
    OutputPointsContainerIterator;

  typedef QuadEdgeMeshPolygonCell<OutputCellType>     OutputPolygonType;
  typedef typename OutputPolygonType::SelfAutoPointer OutputPolygonAutoPointer;
public:
  itkNewMacro( Self );
  itkTypeMacro( QuadEdgeMeshSplitFilter, QuadEdgeMeshToQuadEdgeMeshFilter );

  void SetSeedFaces( const SeedVectorType& iSeeds )
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

  virtual void GenerateData() override;

  OutputPointIdList AddFacePointsToOutputMesh( OutputMeshType* iMesh, InputPolygonType* iPoly );

private:
  QuadEdgeMeshSplitFilter( const Self & );
  void operator =( const Self & );
};
}

#include "itkQuadEdgeMeshSplitFilter.hxx"
#endif
