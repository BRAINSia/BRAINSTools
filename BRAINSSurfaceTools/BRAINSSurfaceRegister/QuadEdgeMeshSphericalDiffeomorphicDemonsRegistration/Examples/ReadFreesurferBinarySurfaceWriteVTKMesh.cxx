/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkFreeSurferBinarySurfaceReaderTest.cxx,v $
  Language:  C++
  Date:      $Date: 2007-02-01 15:22:54 $
  Version:   $Revision: 1.2 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/

#ifdef _MSC_VER
#pragma warning ( disable : 4786 )
#endif

#include "itkQuadEdgeMesh.h"
#include "itkFreeSurferBinarySurfaceReader.h"
#include "itkNormalizeScalarsQuadEdgeMeshFilter.h"
#include "itkQuadEdgeMeshScalarDataVTKPolyDataWriter.h"

#include <iostream>

int main(int argc, char* argv[] )
{
  if( argc < 4 )
    {
    std::cerr << "Usage: " << std::endl;
    std::cerr << argv[0];
    std::cerr << " inputFilename inputDataFilename outputFilename.vtk [normalizeFlag]";
    std::cerr << std::endl;
    return EXIT_FAILURE;
    }

  typedef itk::QuadEdgeMesh<float, 3>                  MeshType;
  typedef itk::FreeSurferBinarySurfaceReader<MeshType> ReaderType;

  ReaderType::Pointer surfaceReader = ReaderType::New();

  surfaceReader->SetFileName( argv[1] );
  surfaceReader->SetDataFileName( argv[2] );

  try
    {
    surfaceReader->Update();
    }
  catch( itk::ExceptionObject & excp )
    {
    std::cerr << "Error during reader Update() " << std::endl;
    std::cerr << excp << std::endl;
    return EXIT_FAILURE;
    }

  MeshType::Pointer mesh = surfaceReader->GetOutput();

  unsigned int numberOfPoints = mesh->GetNumberOfPoints();
  unsigned int numberOfCells  = mesh->GetNumberOfCells();

  std::cout << "numberOfPoints= " << numberOfPoints << std::endl;
  std::cout << "numberOfCells= " << numberOfCells << std::endl;

  typedef itk::NormalizeScalarsQuadEdgeMeshFilter<MeshType> FilterType;

  FilterType::Pointer normalizeFilter = FilterType::New();

  typedef itk::QuadEdgeMeshScalarDataVTKPolyDataWriter<MeshType> WriterType;

  WriterType::Pointer writer = WriterType::New();

  unsigned int normalizationFlag = 0;

  if( argc > 4 )
    {
    normalizationFlag = atoi( argv[4] );
    }

  if( normalizationFlag == 1 )
    {
    normalizeFilter->SetNumberOfIterations( 2 );
    normalizeFilter->SetInput( surfaceReader->GetOutput()  );
    normalizeFilter->Update();
    writer->SetInput( normalizeFilter->GetOutput()  );
    }
  else
    {
    writer->SetInput( surfaceReader->GetOutput()  );
    }

  writer->SetFileName( argv[3] );

  try
    {
    writer->Update();
    }
  catch( itk::ExceptionObject & excp )
    {
    std::cerr << "Error during writer Update() " << std::endl;
    std::cerr << excp << std::endl;
    return EXIT_FAILURE;
    }

  return EXIT_SUCCESS;
}
