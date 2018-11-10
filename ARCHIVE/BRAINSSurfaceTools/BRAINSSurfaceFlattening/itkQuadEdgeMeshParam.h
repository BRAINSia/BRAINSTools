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
#ifndef __itkQuadEdgeMeshParam_h
#define __itkQuadEdgeMeshParam_h

#include "itkQuadEdgeMeshToQuadEdgeMeshFilter.h"
#include "itkQuadEdgeMeshBorderTransform.h"
#include "itkQuadEdgeMeshParamMatrixCoefficients.h"

namespace itk
{
/**
 *  \class QuadEdgeMeshParam
 *
 *  \brief Compute a planar parameterization of the input mesh.
 *  \note Here we have only implemented some parameterizations with fixed
 *        boundary.
 */
template <typename TInputMesh, typename TOutputMesh, typename TSolverTraits>
class QuadEdgeMeshParam :
  public QuadEdgeMeshToQuadEdgeMeshFilter<TInputMesh, TOutputMesh>
{
public:
  ITK_DISALLOW_COPY_AND_ASSIGN(QuadEdgeMeshParam);

  /** Basic types. */
  using Self = QuadEdgeMeshParam;
  using Superclass = QuadEdgeMeshToQuadEdgeMeshFilter<TInputMesh,
                                           TOutputMesh>;
  using Pointer = SmartPointer<Self>;
  using ConstPointer = SmartPointer<const Self>;

  /** Input types. */
  using InputMeshType = TInputMesh;
  using InputMeshPointer = typename InputMeshType::Pointer;
  using InputMeshConstPointer = typename InputMeshType::ConstPointer;
  using InputCoordRepType = typename InputMeshType::CoordRepType;
  using InputPointType = typename InputMeshType::PointType;
  using InputPointVectorType = typename InputPointType::VectorType;
  using InputPointIdentifier = typename InputMeshType::PointIdentifier;
  using InputQEType = typename InputMeshType::QEType;
  using InputVectorType = typename InputMeshType::VectorType;
  using InputEdgeListType = typename InputMeshType::EdgeListType;
  using InputPixelType = typename InputMeshType::PixelType;
  using InputTraits = typename InputMeshType::Traits;

  static constexpr unsigned int InputVDimension = InputMeshType::PointDimension;

  using InputPointsContainer = typename InputMeshType::PointsContainer;
  typedef typename InputMeshType::PointsContainerConstIterator
    InputPointsContainerConstIterator;

  typedef typename InputMeshType::CellsContainerConstIterator
    InputCellsContainerConstIterator;
  using InputEdgeCellType = typename InputMeshType::EdgeCellType;
  using InputPolygonCellType = typename InputMeshType::PolygonCellType;
  using InputPointIdList = typename InputMeshType::PointIdList;

  using InputQEIterator = typename InputQEType::IteratorGeom;

  using InputMapPointIdentifier = std::map<InputPointIdentifier, InputPointIdentifier>;
  typedef typename InputMapPointIdentifier::iterator
    InputMapPoinIdentifierIterator;

  /** Output types. */
  using OutputMeshType = TOutputMesh;
  using OutputMeshPointer = typename OutputMeshType::Pointer;
  using OutputMeshConstPointer = typename OutputMeshType::ConstPointer;
  using OutputCoordRepType = typename OutputMeshType::CoordRepType;
  using OutputPointType = typename OutputMeshType::PointType;
  using OutputPointIdentifier = typename OutputMeshType::PointIdentifier;
  using OutputQEType = typename OutputMeshType::QEType;
  using OutputVectorType = typename OutputMeshType::VectorType;
  using OutputQEIterator = typename OutputQEType::IteratorGeom;
  typedef typename OutputMeshType::PointsContainerIterator
    OutputPointsContainerIterator;

  static constexpr unsigned int OutputVDimension = OutputMeshType::PointDimension;

  using SolverTraits = TSolverTraits;
  using ValueType = typename SolverTraits::ValueType;
  using MatrixType = typename SolverTraits::MatrixType;
  using VectorType = typename SolverTraits::VectorType;

  using MeshBorderTransform = QuadEdgeMeshBorderTransform<
      InputMeshType, InputMeshType>;
  using MeshBorderTransformPointer = typename MeshBorderTransform::Pointer;

  using CoefficientsComputation = MatrixCoefficients<InputMeshType>;
public:

  void SetCoefficientsMethod( CoefficientsComputation* iMethod )
  {
    this->m_CoefficientsMethod = iMethod;
  }

  itkNewMacro( Self );
  itkTypeMacro( MeshParam, QuadEdgeMeshToQuadEdgeMeshFilter );

  itkSetObjectMacro( BorderTransform, MeshBorderTransform );
  itkGetConstObjectMacro( BorderTransform, MeshBorderTransform );
protected:

  QuadEdgeMeshParam();
  virtual ~QuadEdgeMeshParam()
  {
  };

  CoefficientsComputation * m_CoefficientsMethod;

  MeshBorderTransformPointer m_BorderTransform;

  // first is the id of the input mesh and second is the corresponding id
  // in m_BoundaryPtMap
  InputMapPointIdentifier m_BoundaryPtMap;

  // first is the id of the input mesh and second is the corresponding id
  // in m_InternalPtList
  InputMapPointIdentifier m_InternalPtMap;

  std::vector<OutputPointType> m_Border;

  void CopyToOutputBorder();

  /**
   *  \brief From the list of all vertices from the input mesh InputList
   *  and the list of boundary vertices BoundaryList, Store in
   *  m_InternalPtList the list of interior vertices (i.e. vertices in
   *  InputList and not in BoundaryList )
   *
   *  \note I consider ids of points are well chosen (from 0 to
   *        NumberOfPoints)
   */
  void ComputeListOfInteriorVertices();

  /**
   *  \brief Fill matrix iM and vectors Bx and m_By depending on if one
   *  vertex is on the border or not.
   *  \param[in] iM
   *  \param[in,out] ioBx
   *  \param[in,out] ioBy
   */
  void FillMatrix( MatrixType& iM, VectorType& ioBx, VectorType& ioBy );

  /**
   *  \brief Solve linears systems : \f$ iM \cdot oX = iBx \f$ and
   * \f$ iM \cdot oY = iBy \f$
   *
   *  \param[in] iM
   *  \param[in] iBx
   *  \param[in] iBy
   *  \param[out] oX
   *  \param[out] oY
   */
  void SolveLinearSystems( const MatrixType& iM, const VectorType& iBx, const VectorType& iBy, VectorType& oX,
                           VectorType& oY );

  void GenerateData() override;

private:

};
} // end namespace itk

#include "itkQuadEdgeMeshParam.hxx"

#endif
