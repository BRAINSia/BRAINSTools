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
template <typename TInputMesh, typename TOutputMesh>
class QuadEdgeMeshGaussMapFilter :
  public QuadEdgeMeshToQuadEdgeMeshFilter<TInputMesh, TOutputMesh>
{
public:
  typedef QuadEdgeMeshGaussMapFilter Self;
  typedef SmartPointer<Self>         Pointer;
  typedef SmartPointer<const Self>   ConstPointer;
  typedef QuadEdgeMeshToQuadEdgeMeshFilter<TInputMesh, TOutputMesh>
    Superclass;

  /** Run-time type information (and related methods).   */
  itkTypeMacro( QuadEdgeMeshGaussMapFilter, QuadEdgeMeshToQuadEdgeMeshFilter );
  /** New macro for creation of through a Smart Pointer   */
  itkNewMacro( Self );

  /** Input types. */
  typedef TInputMesh                              InputMeshType;
  typedef typename InputMeshType::Pointer         InputMeshPointer;
  typedef typename InputMeshType::ConstPointer    InputMeshConstPointer;
  typedef typename InputMeshType::CoordRepType    InputCoordRepType;
  typedef typename InputMeshType::PointType       InputPointType;
  typedef typename InputPointType::VectorType     InputPointVectorType;
  typedef typename InputMeshType::PointIdentifier InputPointIdentifier;
  typedef typename InputMeshType::QEType          InputQEType;
  typedef typename InputMeshType::VectorType      InputVectorType;
  typedef typename InputMeshType::EdgeListType    InputEdgeListType;
  typedef typename InputMeshType::PixelType       InputPixelType;
  typedef typename InputMeshType::Traits          InputTraits;

  typedef typename InputMeshType::PointsContainer InputPointsContainer;
  typedef typename InputMeshType::PointsContainerConstIterator
    InputPointsContainerConstIterator;

  typedef typename InputMeshType::CellsContainerConstIterator
    InputCellsContainerConstIterator;
  typedef typename InputMeshType::EdgeCellType    InputEdgeCellType;
  typedef typename InputMeshType::PolygonCellType InputPolygonCellType;
  typedef typename InputMeshType::PointIdList     InputPointIdList;

  /** Output types. */
  typedef TOutputMesh                              OutputMeshType;
  typedef typename OutputMeshType::Pointer         OutputMeshPointer;
  typedef typename OutputMeshType::ConstPointer    OutputMeshConstPointer;
  typedef typename OutputMeshType::CoordRepType    OutputCoordRepType;
  typedef typename OutputMeshType::PointType       OutputPointType;
  typedef typename OutputMeshType::PointIdentifier OutputPointIdentifier;
  typedef typename OutputMeshType::QEType          OutputQEType;
  typedef typename OutputMeshType::VectorType      OutputVectorType;
  typedef typename OutputQEType::IteratorGeom      OutputQEIterator;
  typedef typename OutputMeshType::PointsContainerPointer
    OutputPointsContainerPointer;
  typedef typename OutputMeshType::PointsContainerIterator
    OutputPointsContainerIterator;

  static constexpr unsigned int PointDimension = OutputMeshType::PointDimension;

  typedef QuadEdgeMeshExtendedTraits<
      OutputVectorType,
      PointDimension,
      2,
      OutputCoordRepType,
      OutputCoordRepType,
      OutputVectorType,
      bool,
      bool> OutputNormalMeshTraits;
  typedef QuadEdgeMesh<OutputVectorType,
                       PointDimension, OutputNormalMeshTraits> OutputNormalMeshType;
  typedef typename OutputNormalMeshType::Pointer OutputNormalMeshPointer;
  typedef typename OutputNormalMeshType::PointsContainerPointer
    OutputNormalMeshPointsContainerPointer;
  typedef typename OutputNormalMeshType::PointsContainerIterator
    OutputNormalMeshPointsContainerIterator;
  typedef typename OutputNormalMeshType::PointDataContainerPointer
    OutputNormalMeshPointDataContainerPointer;
  typedef typename OutputNormalMeshType::PointDataContainerIterator
    OutputNormalMeshPointDataContainerIterator;

  typedef QuadEdgeMeshNormalFilter<OutputMeshType, OutputNormalMeshType>
    NormalFilterType;
  typedef typename NormalFilterType::Pointer NormalFilterPointer;

  typedef MatrixCoefficients<InputMeshType> CoefficientsComputation;

  void SetCoefficientsMethod( CoefficientsComputation* iMethod )
  {
    (void) iMethod;
  }

  itkSetMacro( Radius, OutputCoordRepType );
protected:
  QuadEdgeMeshGaussMapFilter() : Superclass(), m_Radius(1.)
  {
  }

  ~QuadEdgeMeshGaussMapFilter()
  {
  }

  OutputCoordRepType m_Radius;

  void GenerateData()
  {
    this->CopyInputMeshToOutputMesh();

    OutputMeshPointer output = this->GetOutput();

    NormalFilterPointer normals = NormalFilterType::New();
    normals->SetInput( output );
    normals->SetWeight( NormalFilterType::GOURAUD );
    normals->Update();

    OutputNormalMeshPointer mesh_with_normals = normals->GetOutput();

    OutputNormalMeshPointDataContainerPointer
      pointdata = mesh_with_normals->GetPointData();

    OutputPointIdentifier id( 0 );
    OutputPointType       p;
    unsigned int          dim( 0 );

    OutputNormalMeshPointDataContainerIterator it = pointdata->Begin();
    for( ; it != pointdata->End(); ++it )
      {
      id = it->Index();
      for( dim = 0; dim < PointDimension; ++dim )
        {
        p[dim] = m_Radius * static_cast<OutputCoordRepType>( it.Value()[dim] );
        }
      p.SetEdge( output->FindEdge( id ) );
      output->SetPoint( id, p );
      }
  }

private:
  QuadEdgeMeshGaussMapFilter( const Self & );
  void operator =( const Self & );
};
}
#endif
