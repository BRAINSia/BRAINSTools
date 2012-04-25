/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkRescaleScalarsQuadEdgeMeshFilterTest.cxx,v $
  Language:  C++
  Date:      $Date: 2010-05-17 19:46:31 $
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

#include "itkQuadEdgeMesh.h"
#include "itkRescaleScalarsQuadEdgeMeshFilter.h"

#include "itkQuadEdgeMeshVTKPolyDataReader.h"
#include "itkQuadEdgeMeshScalarDataVTKPolyDataWriter.h"

int main( int argc, char * argv [] )
{
  if( argc < 5 )
    {
    std::cerr << "Missing arguments" << std::endl;
    std::cerr << "Usage: " << std::endl;
    std::cerr << argv[0] << std::endl;
    std::cerr << "inputPolyData outputPolyData";
    std::cerr << "min max";
    std::cerr << std::endl;
    return EXIT_FAILURE;
    }

  typedef double Coord;
  const unsigned int Dimension = 3;

  // Declaration of the type of Mesh
  typedef itk::QuadEdgeMesh<Coord, Dimension> MeshType;

  // Here read a mesh
  typedef itk::QuadEdgeMeshVTKPolyDataReader<MeshType> ReaderType;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[1] );

  typedef itk::RescaleScalarsQuadEdgeMeshFilter<MeshType> RescaleScalarsType;
  RescaleScalarsType::Pointer filter = RescaleScalarsType::New();
  filter->SetInput( reader->GetOutput() );
  filter->SetOutputMinimum(atof(argv[3]) );
  filter->SetOutputMaximum(atof(argv[4]) );
  filter->Update();

  std::cout << "Input min: " << filter->GetInputMinimum() << std::endl;
  std::cout << "Input max: " << filter->GetInputMaximum() << std::endl;

  std::cout << "Output min: " << filter->GetOutputMinimum() << std::endl;
  std::cout << "Output max: " << filter->GetOutputMaximum() << std::endl;

  typedef itk::QuadEdgeMeshScalarDataVTKPolyDataWriter<MeshType> WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetInput( filter->GetOutput() );
  writer->SetFileName(argv[2]);
  writer->Update();

  return EXIT_SUCCESS;
}
