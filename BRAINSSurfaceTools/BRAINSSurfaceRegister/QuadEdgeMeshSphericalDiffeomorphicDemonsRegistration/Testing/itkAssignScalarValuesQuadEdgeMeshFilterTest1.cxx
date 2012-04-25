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

#include "itkAssignScalarValuesQuadEdgeMeshFilter.h"
#include "itkQuadEdgeMeshScalarDataVTKPolyDataWriter.h"
#include "itkQuadEdgeMeshVTKPolyDataReader.h"
#include "itkQuadEdgeMesh.h"

int main( int argc, char * argv [] )
{
  if( argc < 4 )
    {
    std::cerr << "Missing arguments" << std::endl;
    std::cerr << "Usage: " << std::endl;
    std::cerr << argv[0] << std::endl;
    std::cerr << "inputPolyData sourcePolyData outputPolyData";
    std::cerr << std::endl;
    return EXIT_FAILURE;
    }

  typedef float MeshPixelType;
  const unsigned int Dimension = 3;

  typedef itk::QuadEdgeMesh<MeshPixelType, Dimension> InputMeshType;
  typedef itk::QuadEdgeMesh<MeshPixelType, Dimension> SourceMeshType;
  typedef itk::QuadEdgeMesh<MeshPixelType, Dimension> OutputMeshType;

  typedef itk::QuadEdgeMeshVTKPolyDataReader<InputMeshType>  InputReaderType;
  typedef itk::QuadEdgeMeshVTKPolyDataReader<SourceMeshType> SourceReaderType;

  InputReaderType::Pointer inputReader = InputReaderType::New();
  inputReader->SetFileName( argv[1] );

  SourceReaderType::Pointer sourceReader = SourceReaderType::New();
  sourceReader->SetFileName( argv[2] );

  try
    {
    inputReader->Update();
    sourceReader->Update();
    }
  catch( itk::ExceptionObject & exp )
    {
    std::cerr << "Exception thrown while reading the input file " << std::endl;
    std::cerr << exp << std::endl;
    return EXIT_FAILURE;
    }

  typedef itk::AssignScalarValuesQuadEdgeMeshFilter<
      InputMeshType,
      SourceMeshType,
      OutputMeshType>    FilterType;

  FilterType::Pointer filter  = FilterType::New();

  filter->SetInputMesh(inputReader->GetOutput() );
  filter->SetSourceMesh(sourceReader->GetOutput() );
  filter->Update();

  typedef itk::QuadEdgeMeshScalarDataVTKPolyDataWriter<OutputMeshType> WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetInput( filter->GetOutput() );
  writer->SetFileName(argv[3]);
  writer->Update();

  return EXIT_SUCCESS;
}
