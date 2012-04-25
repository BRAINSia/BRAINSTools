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

#include "itkNormalizeScalarsQuadEdgeMeshFilter.h"
#include "itkQuadEdgeMesh.h"
#include "itkQuadEdgeMeshVTKPolyDataReader.h"
#include "itkQuadEdgeMeshScalarDataVTKPolyDataWriter.h"
#include "itkTestingMacros.h"
#include "itkFilterWatcher.h"

int main( int argc, char *argv[] )
{
  if( argc < 3 )
    {
    std::cerr << "Missing Parameters " << std::endl;
    std::cerr << "Usage: " << argv[0];
    std::cerr << " inputMesh outputNormalizedMesh " << std::endl;
    return EXIT_FAILURE;
    }

  typedef float MeshPixelType;
  const unsigned int Dimension = 3;

  typedef itk::QuadEdgeMesh<MeshPixelType, Dimension> MeshType;

  typedef itk::NormalizeScalarsQuadEdgeMeshFilter<MeshType> FilterType;

  FilterType::Pointer normalizeFilter = FilterType::New();

  std::cout << normalizeFilter->GetNameOfClass() << std::endl;
  normalizeFilter->Print( std::cout );

  unsigned int numberOfIterations = 7;

  normalizeFilter->SetNumberOfIterations( numberOfIterations );

  TEST_SET_GET_VALUE( numberOfIterations, normalizeFilter->GetNumberOfIterations() );

  numberOfIterations = 2;

  normalizeFilter->SetNumberOfIterations( numberOfIterations );

  TEST_SET_GET_VALUE( numberOfIterations, normalizeFilter->GetNumberOfIterations() );

  typedef itk::QuadEdgeMeshVTKPolyDataReader<MeshType> ReaderType;

  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[1] );

  FilterWatcher watcher( normalizeFilter, "Normalize Filter");

  typedef itk::QuadEdgeMeshScalarDataVTKPolyDataWriter<MeshType> WriterType;
  WriterType::Pointer writer = WriterType::New();

  writer->SetFileName( argv[2] );

  normalizeFilter->SetInput( reader->GetOutput() );
  writer->SetInput( normalizeFilter->GetOutput() );

  try
    {
    reader->Update();
    normalizeFilter->Update();
    writer->Update();
    }
  catch( itk::ExceptionObject & excp )
    {
    std::cerr << excp << std::endl;
    return EXIT_FAILURE;
    }

  return EXIT_SUCCESS;
}
