/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkMeshToMeshRegistrationMethod.txx,v $
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
#include "itkMeshToMeshRegistrationMethod.h"
#include "itkMeanSquaresMeshToMeshMetric.h"
#include "itkVersorTransform.h"
#include "itkVersorTransformOptimizer.h"
#include "itkNearestNeighborInterpolateMeshFunction.h"
#include "itkTestingMacros.h"

namespace itk
{
template <class TMeshFixed, class TMeshMoving>
class RegistrationMethodHelper : public MeshToMeshRegistrationMethod<TMeshFixed, TMeshMoving>
{
public:
  typedef RegistrationMethodHelper                              Self;
  typedef MeshToMeshRegistrationMethod<TMeshFixed, TMeshMoving> Superclass;
  typedef SmartPointer<Self>                                    Pointer;
  typedef SmartPointer<const Self>                              ConstPointer;

  itkNewMacro( Self );

  itkTypeMacro( RegistrationMethodHelper, MeshToMeshRegistrationMethod );
protected:
  RegistrationMethodHelper()
  {
  }

  virtual ~RegistrationMethodHelper()
  {
  }

private:
};
}

int main( int argc, char * argv [] )
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

  typedef itk::RegistrationMethodHelper<MeshType, MeshType> RegistrationMethodType;

  RegistrationMethodType::Pointer registrator = RegistrationMethodType::New();

  typedef RegistrationMethodType::Superclass RegistrationMethodSuperclassType;

  std::cout << registrator->RegistrationMethodSuperclassType::GetNameOfClass() << std::endl;
  registrator->RegistrationMethodSuperclassType::Print( std::cout );

  std::cout << registrator->GetNameOfClass() << std::endl;
  registrator->Print( std::cout );

  std::cout << std::endl;
  std::cout << "Exercise Exception Catching" << std::endl;
  std::cout << std::endl;

  TRY_EXPECT_EXCEPTION( registrator->Initialize() );

  registrator->SetFixedMesh( meshFixed );

  TRY_EXPECT_EXCEPTION( registrator->Initialize() );

  TEST_SET_GET( meshFixed, registrator->GetFixedMesh() );

  registrator->SetMovingMesh( meshMoving );

  TRY_EXPECT_EXCEPTION( registrator->Initialize() );

  TEST_SET_GET( meshMoving, registrator->GetMovingMesh() );

  typedef itk::MeanSquaresMeshToMeshMetric<MeshType, MeshType> MetricType;
  MetricType::Pointer metric = MetricType::New();

  registrator->SetMetric( metric );

  TRY_EXPECT_EXCEPTION( registrator->Initialize() );

  TEST_SET_GET( metric, registrator->GetMetric() );

  typedef itk::VersorTransformOptimizer OptimizerType;
  OptimizerType::Pointer optimizer = OptimizerType::New();

  registrator->SetOptimizer( optimizer );

  TRY_EXPECT_EXCEPTION( registrator->Initialize() );

  TEST_SET_GET( optimizer, registrator->GetOptimizer() );

  typedef itk::VersorTransform<MetricType::TransformComputationType> TransformType;
  TransformType::Pointer transform = TransformType::New();

  registrator->SetTransform( transform );

  TRY_EXPECT_EXCEPTION( registrator->Initialize() );

  TEST_SET_GET( transform, registrator->GetTransform() );

  typedef itk::NearestNeighborInterpolateMeshFunction<MeshType> InterpolatorType;

  InterpolatorType::Pointer interpolator = InterpolatorType::New();

  registrator->SetInterpolator( interpolator );

  TRY_EXPECT_EXCEPTION( registrator->Initialize() );

  TEST_SET_GET( interpolator, registrator->GetInterpolator() );

  const unsigned int            wrongNumberOfParameters = transform->GetNumberOfParameters() - 1;
  TransformType::ParametersType parameters0( wrongNumberOfParameters );

  registrator->SetInitialTransformParameters( parameters0 );

  TRY_EXPECT_EXCEPTION( registrator->Initialize() );

  TransformType::ParametersType parameters1 = transform->GetParameters();

  registrator->SetInitialTransformParameters( parameters1 );

  TRY_EXPECT_NO_EXCEPTION( registrator->Initialize() );

  TransformType::ParametersType parameters2 = registrator->GetInitialTransformParameters();

  const unsigned int numberOfParameters = transform->GetNumberOfParameters();

  const double tolerance = 1e5;
  for( unsigned int p = 0; p < numberOfParameters; p++ )
    {
    if( vnl_math_abs( parameters2[p] - parameters1[p] ) > tolerance )
      {
      std::cerr << "Error in SetInitialTransformParameters()/GetInitialTransformParameters() " << std::endl;
      return EXIT_FAILURE;
      }
    }

  registrator->SetTransform( NULL );

  TRY_EXPECT_EXCEPTION( registrator->StartRegistration() );

  TRY_EXPECT_EXCEPTION( registrator->MakeOutput(1) );
  TRY_EXPECT_NO_EXCEPTION( registrator->MakeOutput(0) );

  return EXIT_SUCCESS;
}
