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
#ifndef __itkQuadEdgeMeshSphericalParameterization_h
#define __itkQuadEdgeMeshSphericalParameterization_h

#include <itkQuadEdgeMeshToQuadEdgeMeshFilter.h>
#include <itkQuadEdgeMeshParamMatrixCoefficients.h>

namespace itk
{
/**
 * \class QuadEdgeMeshSphericalParameterization
 * \brief Compute the spherical parameterization of the input mesh. The input
 * mesh must be homeomorph to a sphere.
 * \note No test is perform on the topology of the input mesh.
 */
template < typename TInputMesh, typename TOutputMesh, typename TInitializationFilter >
class QuadEdgeMeshSphericalParameterization : public QuadEdgeMeshToQuadEdgeMeshFilter< TInputMesh, TOutputMesh >
{
public:
  using Self = QuadEdgeMeshSphericalParameterization;
  using Superclass = QuadEdgeMeshToQuadEdgeMeshFilter< TInputMesh, TOutputMesh >;
  using Pointer = SmartPointer< Self >;
  using ConstPointer = SmartPointer< const Self >;

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

  typedef typename InputMeshType::PointsContainerPointer       InputPointsContainerPointer;
  typedef typename InputMeshType::PointsContainerConstPointer  InputPointsContainerConstPointer;
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

  using CoefficientsComputation = MatrixCoefficients< InputMeshType >;

  using InitializationFilterType = TInitializationFilter;
  typedef typename InitializationFilterType::Pointer InitializationFilterPointer;
  using PointPairType = std::pair< OutputPointIdentifier, OutputPointIdentifier >;

public:
  void
  SetCoefficientsMethod( CoefficientsComputation * iMethod )
  {
    m_CoefficientsMethod = iMethod;
  }

  itkNewMacro( Self );
  itkTypeMacro( QuadEdgeMeshSphericalParameterization, QuadEdgeMeshToQuadEdgeMeshFilter );

  itkSetMacro( IterationSTOP, unsigned int );
  itkSetMacro( TimeStep, OutputCoordRepType );
  itkSetMacro( Threshold, OutputCoordRepType );
  itkSetMacro( Radius, OutputCoordRepType );

protected:
  QuadEdgeMeshSphericalParameterization()
    : Superclass()
    , m_Iteration( 0 )
    , m_IterationSTOP( 100000 )
    , m_OldEnergy( 1e12 )
    , m_NewEnergy( 1e12 )
    , m_Threshold( 1e-10 )
    , m_TimeStep( 0.1 )
    , m_Radius( 1. )
  {
    m_InitFilter = InitializationFilterType::New();
  }

  ~QuadEdgeMeshSphericalParameterization() {}

  unsigned int       m_Iteration;
  unsigned int       m_IterationSTOP;
  OutputCoordRepType m_OldEnergy;
  OutputCoordRepType m_NewEnergy;

  std::map< PointPairType, OutputCoordRepType > m_CoefficientMap;

  OutputCoordRepType m_Threshold;
  OutputCoordRepType m_TimeStep;
  OutputCoordRepType m_Radius;

  InitializationFilterPointer m_InitFilter;

  CoefficientsComputation * m_CoefficientsMethod;

  void
  ComputeCoefficientMap()
  {
    InputMeshConstPointer input = this->GetInput();

    InputPointsContainerConstPointer points = input->GetPoints();
    InputPointIdentifier             p_id( 0 );
    InputQEType *                    qe( 0 );
    InputQEType *                    qe_it( 0 );
    InputPointType                   p;
    OutputCoordRepType               coeff( 0. );

    for ( InputPointsContainerConstIterator it = points->Begin(); it != points->End(); ++it )
    {
      p_id = it->Index();
      p = it->Value();

      qe = p.GetEdge();
      if ( qe != 0 )
      {
        qe_it = qe;

        do
        {
          if ( p_id < qe_it->GetDestination() )
          {
            coeff = static_cast< OutputCoordRepType >( ( *m_CoefficientsMethod )( input, qe_it ) );
            m_CoefficientMap[PointPairType( p_id, qe_it->GetDestination() )] = coeff;
          }
          qe_it = qe_it->GetOnext();
        } while ( qe_it != qe );
      }
    }
  }

  void
  GenerateData()
  {
    assert( m_CoefficientsMethod != 0 );

    m_Iteration = 0;

    ComputeCoefficientMap();

    m_InitFilter->SetInput( this->GetInput() );
    m_InitFilter->SetRadius( m_Radius );
    m_InitFilter->SetCoefficientsMethod( m_CoefficientsMethod );
    m_InitFilter->Update();

    OutputMeshPointer output = this->GetOutput();
    output->Graft( m_InitFilter->GetOutput() );

    if ( m_IterationSTOP > 0 )
    {
      do
      {
        m_OldEnergy = m_NewEnergy;
        ComputeEnergyAndRelaxVertexLocation();
        // if( m_Iteration % 10 == 0 )
        //           std::cout <<m_Iteration <<" " <<m_NewEnergy <<" "
        //             <<itk::Math::abs ( m_OldEnergy - m_NewEnergy )
        //             <<std::endl;
        ++m_Iteration;
      } while ( itk::Math::abs( m_OldEnergy - m_NewEnergy ) > m_Threshold && m_Iteration < m_IterationSTOP );

      std::cout << m_Iteration << std::endl;
    }
  }

  inline OutputCoordRepType
  ComputeQuadEdgeEnergy( const OutputCoordRepType & iCoeff, const OutputPointType & iOrg,
                         const OutputPointType & iDest )
  {
    return iCoeff * iOrg.SquaredEuclideanDistanceTo( iDest );
  }

  inline OutputVectorType
  ComputeQuadEdgeAbsoluteDerivative( const OutputCoordRepType & iCoeff, const OutputPointType & iOrg,
                                     const OutputPointType & iDest )
  {
    return iCoeff * ( iDest - iOrg );
  }

  /** E = \sum_K k_{uv} \| n_u - n_v \|^2 */
  OutputVectorType
  ComputeEnergyAndAbsoluteDerivative( const OutputPointIdentifier & iId )
  {
    OutputMeshPointer output = this->GetOutput();
    OutputQEType *    qe = output->FindEdge( iId );
    OutputVectorType  oVector;

    oVector.Fill( 0. );
    OutputPointIdentifier id_dest;

    if ( qe != 0 )
    {
      OutputQEType * qe_it = qe;

      OutputPointType p_org = output->GetPoint( iId );
      OutputPointType p_dest;

      OutputVectorType delta;
      delta.Fill( 0. );

      OutputVectorType   n = p_org.GetVectorFromOrigin();
      OutputCoordRepType coeff( 0. );

      do
      {
        id_dest = qe_it->GetDestination();
        p_dest = output->GetPoint( id_dest );

        coeff = ( iId < id_dest ) ? m_CoefficientMap[PointPairType( iId, id_dest )]
                                  : m_CoefficientMap[PointPairType( id_dest, iId )];

        m_NewEnergy += ComputeQuadEdgeEnergy( coeff, p_org, p_dest );
        delta += ComputeQuadEdgeAbsoluteDerivative( coeff, p_org, p_dest );
        qe_it = qe_it->GetOnext();
      } while ( qe_it != qe );

      oVector = delta - ( delta * n ) * n;
    }

    return oVector;
  }

  void
  ComputeEnergyAndRelaxVertexLocation()
  {
    m_NewEnergy = static_cast< OutputCoordRepType >( 0. );
    OutputMeshPointer output = this->GetOutput();

    OutputPointsContainerPointer points = output->GetPoints();
    OutputVectorType             v;
    OutputPointType              p, q;
    OutputPointIdentifier        id( 0 );
    unsigned int                 dim( 0 );
    OutputCoordRepType           norm( 0. );
    for ( OutputPointsContainerIterator it = points->Begin(); it != points->End(); ++it )
    {
      id = it->Index();
      p = it->Value();

      if ( p.GetEdge() != 0 )
      {
        v = ComputeEnergyAndAbsoluteDerivative( id );
        norm = 0.;
        for ( dim = 0; dim < PointDimension; ++dim )
        {
          p[dim] = ( 1. - m_TimeStep ) * p[dim] + m_TimeStep * v[dim];
          norm += p[dim] * p[dim];
        }

        norm = m_Radius / std::sqrt( norm );
        for ( dim = 0; dim < PointDimension; ++dim )
        {
          p[dim] *= norm;
        }

        output->SetPoint( id, p );
        q = output->GetPoint( id );
      }
      else
      {
        std::cout << "p.GetEdge()==0" << std::endl;
      }
    }
  }

private:
  QuadEdgeMeshSphericalParameterization( const Self & );
  void
  operator=( const Self & );
};
} // namespace itk
#endif
