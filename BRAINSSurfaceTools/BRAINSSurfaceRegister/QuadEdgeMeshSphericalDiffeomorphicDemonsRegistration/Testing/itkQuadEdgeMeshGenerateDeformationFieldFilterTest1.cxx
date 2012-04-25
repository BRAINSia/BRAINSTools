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

#include "itkQuadEdgeMeshGenerateDeformationFieldFilter.h"
#include "itkQuadEdgeMesh.h"
#include "itkQuadEdgeMeshVTKPolyDataReader.h"
#include "itkQuadEdgeMeshVectorDataVTKPolyDataWriter.h"
#include "itkFilterWatcher.h"

int main( int argc, char *argv[] )
{
  if( argc < 4 )
    {
    std::cerr << "Missing Parameters " << std::endl;
    std::cerr << "Usage: " << argv[0];
    std::cerr << " fixedMeshFile  deformedMeshFile  outputMeshDisplacementFieldfile " << std::endl;
    return EXIT_FAILURE;
    }

  typedef float MeshPixelType;
  const unsigned int Dimension = 3;

  typedef itk::QuadEdgeMesh<MeshPixelType, Dimension> FixedMeshType;
  typedef itk::QuadEdgeMesh<MeshPixelType, Dimension> RegisteredMeshType;

  typedef itk::Vector<float, Dimension> DeformationVectorType;

  typedef itk::QuadEdgeMeshTraits<DeformationVectorType, Dimension, bool, bool> VectorPointSetTraits;

  typedef itk::QuadEdgeMesh<
      DeformationVectorType, Dimension, VectorPointSetTraits>   DeformationFieldMeshType;

  typedef itk::QuadEdgeMeshGenerateDeformationFieldFilter<
      FixedMeshType, RegisteredMeshType, DeformationFieldMeshType>   DeformationFilterType;

  DeformationFilterType::Pointer deformationFilter = DeformationFilterType::New();

  std::cout << deformationFilter->GetNameOfClass() << std::endl;
  deformationFilter->Print( std::cout );

  typedef itk::QuadEdgeMeshVTKPolyDataReader<FixedMeshType>      FixedReaderType;
  typedef itk::QuadEdgeMeshVTKPolyDataReader<RegisteredMeshType> RegisteredReaderType;

  FixedReaderType::Pointer fixedReader = FixedReaderType::New();
  fixedReader->SetFileName( argv[1] );

  RegisteredReaderType::Pointer registeredReader = RegisteredReaderType::New();
  registeredReader->SetFileName( argv[2] );

  try
    {
    fixedReader->Update();
    registeredReader->Update();
    }
  catch( itk::ExceptionObject & exp )
    {
    std::cerr << "Exception thrown while reading the input file " << std::endl;
    std::cerr << exp << std::endl;
    return EXIT_FAILURE;
    }

  deformationFilter->SetInputMesh( fixedReader->GetOutput() );
  deformationFilter->SetDestinationPoints( registeredReader->GetOutput() );

  FilterWatcher watcher( deformationFilter, "Demons Filter");

  try
    {
    deformationFilter->Update();
    }
  catch( itk::ExceptionObject & exp )
    {
    std::cerr << "Exception thrown while running the Demons filter " << std::endl;
    std::cerr << exp << std::endl;
    return EXIT_FAILURE;
    }

  typedef itk::QuadEdgeMeshVectorDataVTKPolyDataWriter<DeformationFieldMeshType> WriterType;
  WriterType::Pointer writer = WriterType::New();

  writer->SetFileName( argv[3] );
  writer->SetInput( deformationFilter->GetOutput() );

  try
    {
    writer->Update();
    }
  catch( itk::ExceptionObject & excp )
    {
    std::cerr << excp << std::endl;
    return EXIT_FAILURE;
    }

  std::cout << std::endl;
  std::cout << "Testing second call to Update(), the filter shouldn't run again" << std::endl;
  deformationFilter->Update();

  return EXIT_SUCCESS;
}
