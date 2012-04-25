/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkMeshToMeshMetric.txx,v $
  Language:  C++
  Date:      $Date: 2003-11-08 17:58:32 $
  Version:   $Revision: 1.1 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/

#include "itkQuadEdgeMesh.h"
#include "itkVTKPolyDataReader.h"
#include "itkQuadEdgeMeshScalarDataVTKPolyDataWriter.h"
#include "itkMeshToMeshMetric.h"
#include "itkNearestNeighborInterpolateMeshFunction.h"
#include "itkVersorTransform.h"
#include "itkTestingMacros.h"

namespace itk
{
template <class TMeshFixed, class TMeshMoving>
class MetricHelper : public MeshToMeshMetric<TMeshFixed, TMeshMoving>
{
public:
  typedef MetricHelper                              Self;
  typedef MeshToMeshMetric<TMeshFixed, TMeshMoving> Superclass;
  typedef SmartPointer<Self>                        Pointer;
  typedef SmartPointer<const Self>                  ConstPointer;

  itkNewMacro( Self );

  itkTypeMacro( MetricHelper, MeshToMeshMetric );

  typedef typename Superclass::MeasureType             MeasureType;
  typedef typename Superclass::DerivativeType          DerivativeType;
  typedef typename Superclass::TransformParametersType TransformParametersType;
protected:
  MetricHelper()
  {
  }

  virtual ~MetricHelper()
  {
  }

  /**  Get the value for single valued optimizers. */
  MeasureType GetValue( const TransformParametersType & itkNotUsed(parameters) ) const
  {
    return 1.0;
  }

  /**  Get derivatives for multiple valued optimizers. */
  void GetDerivative( const TransformParametersType & itkNotUsed(parameters),
                      DerivativeType& derivative ) const
  {
    derivative.Fill( 0.0 );
  }

  /**  Get value and derivatives for multiple valued optimizers. */
  void GetValueAndDerivative( const TransformParametersType & itkNotUsed(parameters),
                              MeasureType& Value, DerivativeType& derivative ) const
  {
    Value = 1.0;
    derivative.Fill( 0.0 );
  }
};
}

int main( int argc, char* * argv )
{
  if( argc != 2 )
    {
    std::cout << "It requires 1 arguments" << std::endl;
    std::cout << "1-Input FileName" << std::endl;
    return EXIT_FAILURE;
    }

  typedef double Coord;

  typedef itk::QuadEdgeMesh<Coord, 3>      MeshType;
  typedef MeshType::Pointer                MeshPointer;
  typedef itk::VTKPolyDataReader<MeshType> ReaderType;

  typedef itk::QuadEdgeMeshScalarDataVTKPolyDataWriter<MeshType> WriterType;

  ReaderType::Pointer readerFixed = ReaderType::New();
  readerFixed->SetFileName( argv[1] );

  ReaderType::Pointer readerMoving = ReaderType::New();
  readerMoving->SetFileName( argv[1] );

  try
    {
    readerFixed->Update();
    readerMoving->Update();
    }
  catch( itk::ExceptionObject & exp )
    {
    std::cerr << "Exception thrown while reading the input file " << std::endl;
    std::cerr << exp << std::endl;
    return EXIT_FAILURE;
    }

  MeshPointer meshFixed  = readerFixed->GetOutput();
  MeshPointer meshMoving = readerMoving->GetOutput();

  typedef itk::MetricHelper<MeshType, MeshType> MetricType;
  typedef MetricType::Superclass                MetricSuperclassType;

  MetricType::Pointer metric = MetricType::New();

  TRY_EXPECT_EXCEPTION( metric->Initialize() );

  typedef itk::VersorTransform<MetricType::TransformComputationType> TransformType;
  TransformType::Pointer transform = TransformType::New();

  metric->SetTransform( transform );

  TRY_EXPECT_EXCEPTION( metric->Initialize() );

  typedef itk::NearestNeighborInterpolateMeshFunction<MeshType> InterpolatorType;

  InterpolatorType::Pointer interpolator = InterpolatorType::New();

  metric->SetInterpolator( interpolator );

  TRY_EXPECT_EXCEPTION( metric->Initialize() );

  MetricType::InterpolatorType::ConstPointer interpolator2 =
    metric->GetInterpolator();

  if( interpolator2.GetPointer() != interpolator.GetPointer() )
    {
    std::cerr << "Error in Set/GetInterpolator() " << std::endl;
    return EXIT_FAILURE;
    }

  metric->SetMovingMesh( meshMoving );

  TRY_EXPECT_EXCEPTION( metric->Initialize() );

  metric->SetFixedMesh( meshFixed );

  TRY_EXPECT_NO_EXCEPTION( metric->Initialize() );

  metric->SetTransform( NULL );

  TRY_EXPECT_EXCEPTION( metric->SetTransformParameters( transform->GetParameters() ) );

  metric->SetTransform( transform );

  TRY_EXPECT_NO_EXCEPTION( metric->SetTransformParameters( transform->GetParameters() ) );

  std::cout << metric->MetricSuperclassType::GetNameOfClass() << std::endl;
  metric->MetricSuperclassType::Print( std::cout );

  std::cout << metric->GetNameOfClass() << std::endl;
  metric->Print( std::cout );

  // Exercising error management of metric verification for the presence of
  // point coordinates and point data.
  MeshType::PointDataContainer::Pointer pointData = meshFixed->GetPointData();
  meshFixed->SetPointData( NULL );
  TRY_EXPECT_EXCEPTION( metric->Initialize() );
  meshFixed->SetPointData( pointData );
  TRY_EXPECT_NO_EXCEPTION( metric->Initialize() );

  MeshType::PointsContainer::Pointer points = meshFixed->GetPoints();
  meshFixed->SetPoints( NULL );
  TRY_EXPECT_EXCEPTION( metric->Initialize() );
  meshFixed->SetPoints( points );
  TRY_EXPECT_NO_EXCEPTION( metric->Initialize() );

  pointData = meshMoving->GetPointData();
  meshMoving->SetPointData( NULL );
  TRY_EXPECT_EXCEPTION( metric->Initialize() );
  meshMoving->SetPointData( pointData );
  TRY_EXPECT_NO_EXCEPTION( metric->Initialize() );

  points = meshMoving->GetPoints();
  meshMoving->SetPoints( NULL );
  TRY_EXPECT_EXCEPTION( metric->Initialize() );
  meshMoving->SetPoints( points );
  TRY_EXPECT_NO_EXCEPTION( metric->Initialize() );

  return EXIT_SUCCESS;
}
