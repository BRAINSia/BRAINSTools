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
#ifndef __itkQuadEdgeMeshRayTracingFilter_h
#define __itkQuadEdgeMeshRayTracingFilter_h

#include <itkQuadEdgeMeshToQuadEdgeMeshFilter.h>
#include <itkQuadEdgeMeshParamMatrixCoefficients.h>

namespace itk
{
/**
 * \class QuadEdgeMeshRayTracingFilter
 * \brief This filter first computes the center of mass C of the input mesh.
 * Then it projects each vertex on a sphere centered at C.
 */
template <typename TInputMesh, typename TOutputMesh>
class QuadEdgeMeshRayTracingFilter : public QuadEdgeMeshToQuadEdgeMeshFilter<TInputMesh, TOutputMesh>
{
public:
  using Self = QuadEdgeMeshRayTracingFilter;
  using Pointer = SmartPointer<Self>;
  using ConstPointer = SmartPointer<const Self>;
  using Superclass = QuadEdgeMeshToQuadEdgeMeshFilter<TInputMesh, TOutputMesh>;

  /** Run-time type information (and related methods).   */
  itkTypeMacro(QuadEdgeMeshRayTracingFilter, QuadEdgeMeshToQuadEdgeMeshFilter);
  /** New macro for creation of through a Smart Pointer   */
  itkNewMacro(Self);

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

  using InputPointsContainer = typename InputMeshType::PointsContainer;
  typedef typename InputMeshType::PointsContainerConstIterator InputPointsContainerConstIterator;

  typedef typename InputMeshType::CellsContainerConstIterator InputCellsContainerConstIterator;
  using InputEdgeCellType = typename InputMeshType::EdgeCellType;
  using InputPolygonCellType = typename InputMeshType::PolygonCellType;
  using InputPointIdList = typename InputMeshType::PointIdList;

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
  typedef typename OutputMeshType::PointsContainerPointer  OutputPointsContainerPointer;
  typedef typename OutputMeshType::PointsContainerIterator OutputPointsContainerIterator;

  static constexpr unsigned int PointDimension = OutputMeshType::PointDimension;

  using CoefficientsComputation = MatrixCoefficients<InputMeshType>;

  void
  SetCoefficientsMethod(CoefficientsComputation * iMethod)
  {
    (void)iMethod;
  }

  itkSetMacro(Radius, OutputCoordRepType);

protected:
  QuadEdgeMeshRayTracingFilter()
    : Superclass()
    , m_Radius(1.)
  {}

  ~QuadEdgeMeshRayTracingFilter() {}

  OutputPointType    m_Center;
  OutputCoordRepType m_Radius;

  void
  GenerateData()
  {
    this->CopyInputMeshToOutputMesh();

    OutputMeshPointer             output = this->GetOutput();
    OutputPointsContainerPointer  points = output->GetPoints();
    OutputPointsContainerIterator p_it = points->Begin();

    OutputCoordRepType norm2;
    OutputPointType    u;
    unsigned int       dim;
    for (p_it = points->Begin(); p_it != points->End(); ++p_it)
    {
      norm2 = 0.;
      u.SetEdge(p_it->Value().GetEdge());
      for (dim = 0; dim < PointDimension; ++dim)
      {
        u[dim] = p_it->Value()[dim] - m_Center[dim];
        norm2 += u[dim] * u[dim];
      }
      norm2 = m_Radius / std::sqrt(norm2);
      for (dim = 0; dim < PointDimension; ++dim)
      {
        u[dim] *= norm2;
      }
      points->SetElement(p_it->Index(), u);
    }
  }

  /** \brief Compute the Center of mass of the input mesh. */
  void
  ComputeCenterOfMass()
  {
    m_Center.Fill(0.);

    OutputMeshPointer             output = this->GetOutput();
    OutputPointsContainerPointer  points = output->GetPoints();
    OutputPointsContainerIterator p_it = points->Begin();

    unsigned int  dim;
    unsigned long k = 0;
    for (; p_it != points->End(); ++p_it, ++k)
    {
      for (dim = 0; dim < PointDimension; ++dim)
      {
        m_Center[dim] += p_it->Value()[dim];
      }
    }

    OutputCoordRepType inv = 1. / static_cast<OutputCoordRepType>(k);
    for (dim = 0; dim < PointDimension; ++dim)
    {
      m_Center[dim] *= inv;
    }
  }

private:
  QuadEdgeMeshRayTracingFilter(const Self &);
  void
  operator=(const Self &);
};
} // namespace itk
#endif
