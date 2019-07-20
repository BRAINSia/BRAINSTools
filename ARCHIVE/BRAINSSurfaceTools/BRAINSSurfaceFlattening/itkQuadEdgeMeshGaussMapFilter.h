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
#ifndef __itkQuadEdgeMeshGaussMapFilter_h
#define __itkQuadEdgeMeshGaussMapFilter_h

#include <itkQuadEdgeMeshExtendedTraits.h>
#include <itkQuadEdgeMeshNormalFilter.h>
#include <itkQuadEdgeMeshToQuadEdgeMeshFilter.h>
#include <itkQuadEdgeMeshParamMatrixCoefficients.h>

namespace itk
{
/**
 * \class QuadEdgeMeshGaussMapFilter
 * \brief Compute the Gauss map of a given mesh. It change the position of
 * each vertex to its computed normal.
 */
template < typename TInputMesh, typename TOutputMesh >
class QuadEdgeMeshGaussMapFilter : public QuadEdgeMeshToQuadEdgeMeshFilter< TInputMesh, TOutputMesh >
{
public:
  using Self = QuadEdgeMeshGaussMapFilter;
  using Pointer = SmartPointer< Self >;
  using ConstPointer = SmartPointer< const Self >;
  using Superclass = QuadEdgeMeshToQuadEdgeMeshFilter< TInputMesh, TOutputMesh >;

  /** Run-time type information (and related methods).   */
  itkTypeMacro( QuadEdgeMeshGaussMapFilter, QuadEdgeMeshToQuadEdgeMeshFilter );
  /** New macro for creation of through a Smart Pointer   */
  itkNewMacro( Self );

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

  typedef QuadEdgeMeshExtendedTraits< OutputVectorType, PointDimension, 2, OutputCoordRepType, OutputCoordRepType,
                                      OutputVectorType, bool, bool >
    OutputNormalMeshTraits;
  using OutputNormalMeshType = QuadEdgeMesh< OutputVectorType, PointDimension, OutputNormalMeshTraits >;
  using OutputNormalMeshPointer = typename OutputNormalMeshType::Pointer;
  typedef typename OutputNormalMeshType::PointsContainerPointer     OutputNormalMeshPointsContainerPointer;
  typedef typename OutputNormalMeshType::PointsContainerIterator    OutputNormalMeshPointsContainerIterator;
  typedef typename OutputNormalMeshType::PointDataContainerPointer  OutputNormalMeshPointDataContainerPointer;
  typedef typename OutputNormalMeshType::PointDataContainerIterator OutputNormalMeshPointDataContainerIterator;

  using NormalFilterType = QuadEdgeMeshNormalFilter< OutputMeshType, OutputNormalMeshType >;
  using NormalFilterPointer = typename NormalFilterType::Pointer;

  using CoefficientsComputation = MatrixCoefficients< InputMeshType >;

  void
  SetCoefficientsMethod( CoefficientsComputation * iMethod )
  {
    (void)iMethod;
  }

  itkSetMacro( Radius, OutputCoordRepType );

protected:
  QuadEdgeMeshGaussMapFilter()
    : Superclass()
    , m_Radius( 1. )
  {}

  ~QuadEdgeMeshGaussMapFilter() {}

  OutputCoordRepType m_Radius;

  void
  GenerateData()
  {
    this->CopyInputMeshToOutputMesh();

    OutputMeshPointer output = this->GetOutput();

    NormalFilterPointer normals = NormalFilterType::New();
    normals->SetInput( output );
    normals->SetWeight( NormalFilterType::GOURAUD );
    normals->Update();

    OutputNormalMeshPointer mesh_with_normals = normals->GetOutput();

    OutputNormalMeshPointDataContainerPointer pointdata = mesh_with_normals->GetPointData();

    OutputPointIdentifier id( 0 );
    OutputPointType       p;
    unsigned int          dim( 0 );

    OutputNormalMeshPointDataContainerIterator it = pointdata->Begin();
    for ( ; it != pointdata->End(); ++it )
    {
      id = it->Index();
      for ( dim = 0; dim < PointDimension; ++dim )
      {
        p[dim] = m_Radius * static_cast< OutputCoordRepType >( it.Value()[dim] );
      }
      p.SetEdge( output->FindEdge( id ) );
      output->SetPoint( id, p );
    }
  }

private:
  QuadEdgeMeshGaussMapFilter( const Self & );
  void
  operator=( const Self & );
};
} // namespace itk
#endif
