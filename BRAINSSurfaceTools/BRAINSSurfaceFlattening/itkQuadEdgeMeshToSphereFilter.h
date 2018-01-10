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
#ifndef __itkQuadEdgeMeshToSphereFilter_h
#define __itkQuadEdgeMeshToSphereFilter_h

#include <itkQuadEdgeMeshToQuadEdgeMeshFilter.h>
#include "itkQuadEdgeMeshSplitFilter.h"
#include "itkQuadEdgeMeshBoundarySmoothFilter.h"
#include "itkQuadEdgeMeshBorderTransform.h"
#include <itkQuadEdgeMeshParamMatrixCoefficients.h>
#include "itkQuadEdgeMeshParam.h"

namespace itk
{
/**
 * \class QuadEdgeMeshToSphereFilter
 * \brief Map the input mesh to a sphere. First the input mesh is split into
 * two meshes. These two meshes are then parameterized onto a planar disk.
 * Then each disk is mapped to one sphere hemisphere by inverse stereo
 * projection.
*/
template <class TInputMesh,
          class TOutputMesh,
          class TSolverTraits>
class QuadEdgeMeshToSphereFilter :
  public QuadEdgeMeshToQuadEdgeMeshFilter<TInputMesh, TOutputMesh>
{
public:
  typedef QuadEdgeMeshToSphereFilter                                Self;
  typedef SmartPointer<Self>                                        Pointer;
  typedef SmartPointer<const Self>                                  ConstPointer;
  typedef QuadEdgeMeshToQuadEdgeMeshFilter<TInputMesh, TOutputMesh> Superclass;

  /** Run-time type information (and related methods).   */
  itkTypeMacro( QuadEdgeMeshToSphereFilter, QuadEdgeMeshToQuadEdgeMeshFilter );
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
  typedef typename InputMeshType::CellIdentifier  InputCellIdentifier;

  typedef std::vector<InputCellIdentifier> SeedVectorType;

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

  itkStaticConstMacro( PointDimension, unsigned int,
                       OutputMeshType::PointDimension );

  typedef TSolverTraits SolverTraits;

  typedef QuadEdgeMeshSplitFilter<InputMeshType,
                                  OutputMeshType>                                SplitFilterType;
  typedef typename SplitFilterType::Pointer SplitFilterPointer;

  typedef QuadEdgeMeshBoundarySmoothFilter<InputMeshType,
                                           OutputMeshType>
    BoundarySmoothFilterType;
  typedef typename BoundarySmoothFilterType::Pointer BoundarySmoothFilterPointer;

  typedef QuadEdgeMeshBorderTransform<OutputMeshType, OutputMeshType>
    BorderTransformType;
  typedef typename BorderTransformType::Pointer
    BorderTransformPointer;

  typedef QuadEdgeMeshParam<OutputMeshType, OutputMeshType, SolverTraits>
    ParametrizationType;
  typedef typename ParametrizationType::Pointer
    ParametrizationPointer;

  typedef MatrixCoefficients<OutputMeshType> CoefficientsComputation;

  /**
   * \brief Provide matrix coefficients method for planar parameterization.
  */
  void SetCoefficientsMethod( CoefficientsComputation* iMethod )
  {
    this->m_CoefficientsMethod = iMethod;
  }

  void SetSeedFaces( const SeedVectorType& iSeeds )
  {
    m_SeedFaces = iSeeds;
  }

  itkSetMacro( Radius, OutputCoordRepType );
protected:
  QuadEdgeMeshToSphereFilter() : Superclass(),
    m_CoefficientsMethod( nullptr ), m_Radius( 1. )
  {
  }

  ~QuadEdgeMeshToSphereFilter()
  {
  }

  CoefficientsComputation * m_CoefficientsMethod;
  OutputCoordRepType        m_Radius;
  SeedVectorType            m_SeedFaces;

  void GenerateData() ITK_OVERRIDE
  {
    assert( m_CoefficientsMethod != 0 );
    this->CopyInputMeshToOutputMesh();
    OutputMeshPointer output = this->GetOutput();

    // split the input mesh into two meshes
    SplitFilterPointer split_filter = SplitFilterType::New();
    split_filter->SetInput( this->GetInput() );
    split_filter->SetSeedFaces(m_SeedFaces);
    split_filter->Update();

    std::cout << "Split DONE!" << std::endl;

    BoundarySmoothFilterPointer boundary_smooth = BoundarySmoothFilterType::New();
    boundary_smooth->SetInputMesh1(split_filter->GetOutput(0) );
    boundary_smooth->SetInputMesh2(split_filter->GetOutput(1) );
    boundary_smooth->SetIterations(5);
    boundary_smooth->Update();

    BorderTransformPointer border_transform = BorderTransformType::New();
    border_transform->SetInput( boundary_smooth->GetOutputMesh1() );
    border_transform->SetTransformType(
      BorderTransformType::DISK_BORDER_TRANSFORM );
    border_transform->SetRadius( 1. );

    // Inverse stereo projection
    OutputPointsContainerPointer  points;
    OutputPointType               p, q;
    OutputCoordRepType            den, r2;
    OutputPointsContainerIterator p_it;

      {
      ParametrizationPointer param0 = ParametrizationType::New();
      param0->SetInput( boundary_smooth->GetOutputMesh1() );
      param0->SetBorderTransform( border_transform );
      param0->SetCoefficientsMethod( m_CoefficientsMethod );
      param0->Update();

      std::cout << "First part parameterization on a disk: DONE!" << std::endl;

      points = param0->GetOutput()->GetPoints();
      for( p_it = points->Begin(); p_it != points->End(); ++p_it )
        {
        q = output->GetPoint( p_it->Index() );
        p = p_it->Value();
        r2 = p[0] * p[0] + p[1] * p[1];
        den = 1. / ( 1. + r2 );
        q[0] = m_Radius * ( 2. * p[0] * den );
        q[1] = m_Radius * ( 2. * p[1] * den );
        q[2] = m_Radius * ( 2. * r2 * den  - 1. );

        output->SetPoint( p_it->Index(), q );
        }

      std::cout << "Inverse stereo projection Part 1: DONE!" << std::endl;
      }

      {
      ParametrizationPointer param1 = ParametrizationType::New();
      param1->SetInput(boundary_smooth->GetOutputMesh2() );
      param1->SetBorderTransform( border_transform );
      param1->SetCoefficientsMethod( m_CoefficientsMethod );
      param1->Update();

      std::cout << "Second part parameterization on a disk: DONE!" << std::endl;

      points = param1->GetOutput()->GetPoints();
      for( p_it = points->Begin(); p_it != points->End(); ++p_it )
        {
        q = output->GetPoint( p_it->Index() );
        p = p_it->Value();
        r2 = p[0] * p[0] + p[1] * p[1];
        den = 1. / ( 1. + r2 );
        q[0] = m_Radius * ( 2. * p[0] * den );
        q[1] = m_Radius * ( 2. * p[1] * den );
        q[2] = m_Radius * ( -2. * r2 * den + 1. );

        output->SetPoint( p_it->Index(), q );
        }

      std::cout << "Inverse stereo projection Part 2: DONE!" << std::endl;
      }
  }

private:
  QuadEdgeMeshToSphereFilter( const Self & );
  void operator =( const Self & );
};
}
#endif
