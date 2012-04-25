/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: ImageRegistration8.cxx,v $
  Language:  C++
  Date:      $Date: 2008-03-10 19:46:31 $
  Version:   $Revision: 1.37 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif

#include "itkQuadEdgeMeshVectorPixelValuesSmoothingFilter.h"
#include "itkQuadEdgeMeshVectorDataVTKPolyDataWriter.h"
#include "itkQuadEdgeMeshVTKPolyDataReader.h"
#include "itkQuadEdgeMesh.h"

int main( int argc, char *argv[] )
{
  if( argc < 5 )
    {
    std::cerr << "Missing Parameters " << std::endl;
    std::cerr << "Usage: " << argv[0];
    std::cerr << " inputMeshFile  outputMeshfile lambda iterations" << std::endl;
    return EXIT_FAILURE;
    }

  const unsigned int Dimension = 3;
  typedef itk::Vector<float, Dimension> MeshPixelType;

  typedef itk::QuadEdgeMesh<MeshPixelType, Dimension> InputMeshType;
  typedef itk::QuadEdgeMesh<MeshPixelType, Dimension> OutputMeshType;

  typedef itk::QuadEdgeMeshVectorPixelValuesSmoothingFilter<
      InputMeshType, OutputMeshType>                       FilterType;

  FilterType::Pointer filter = FilterType::New();

  typedef itk::QuadEdgeMeshVTKPolyDataReader<InputMeshType> ReaderType;

  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[1] );

  FilterType::OutputPointType sphereCenter;
  sphereCenter.Fill( 0.0 );

  const double sphereRadius = 100.0;

  filter->SetSphereCenter( sphereCenter );
  filter->SetSphereRadius( sphereRadius );

  const double       lambdaValue = atof( argv[3] );
  const unsigned int numberOfIterations = atoi( argv[4] );

  filter->SetLambda( lambdaValue );
  filter->SetMaximumNumberOfIterations( numberOfIterations );

  typedef itk::QuadEdgeMeshVectorDataVTKPolyDataWriter<OutputMeshType> WriterType;

  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( argv[2] );

  filter->SetInput( reader->GetOutput() );
  writer->SetInput( filter->GetOutput() );

  filter->Update();

  try
    {
    writer->Update();
    }
  catch( itk::ExceptionObject & excp )
    {
    std::cerr << excp << std::endl;
    return EXIT_FAILURE;
    }

  return EXIT_SUCCESS;
}
