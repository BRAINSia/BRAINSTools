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

#include "itkQuadEdgeMeshScalarPixelValuesSmoothingFilter.h"
#include "itkQuadEdgeMeshVTKPolyDataReader.h"
#include "itkQuadEdgeMeshScalarDataVTKPolyDataWriter.h"
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
  typedef float MeshPixelType;

  typedef itk::QuadEdgeMesh<MeshPixelType, Dimension> InputMeshType;
  typedef itk::QuadEdgeMesh<MeshPixelType, Dimension> OutputMeshType;

  typedef itk::QuadEdgeMeshVTKPolyDataReader<InputMeshType> ReaderType;

  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[1] );

  typedef itk::QuadEdgeMeshScalarPixelValuesSmoothingFilter<
      InputMeshType, OutputMeshType>                       FilterType;

  FilterType::Pointer filter = FilterType::New();

  const double       lambda = atof( argv[3] );
  const unsigned int numberOfIterations = atoi( argv[4] );

  filter->SetLambda( lambda );
  filter->SetMaximumNumberOfIterations( numberOfIterations );

  typedef itk::QuadEdgeMeshScalarDataVTKPolyDataWriter<OutputMeshType> WriterType;
  WriterType::Pointer writer = WriterType::New();

  filter->SetInput( reader->GetOutput() );
  writer->SetInput( filter->GetOutput() );

  filter->Update();

  writer->SetFileName( argv[2] );

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
