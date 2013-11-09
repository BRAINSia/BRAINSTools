/*=========================================================================

 Program:   BRAINS (Brain Research: Analysis of Images, Networks, and Systems)
 Module:    $RCSfile: $
 Language:  C++
 Date:      $Date: 2011/07/09 14:53:40 $
 Version:   $Revision: 1.0 $

 Copyright (c) University of Iowa Department of Radiology. All rights reserved.
 See GTRACT-Copyright.txt or http://mri.radiology.uiowa.edu/copyright/GTRACT-Copyright.txt
 for details.

 This software is distributed WITHOUT ANY WARRANTY; without even
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE.  See the above copyright notices for more information.

 =========================================================================*/

#include "itkQuadEdgeMeshScalarPixelValuesSmoothingFilter.h"
#include "itkQuadEdgeMeshVTKPolyDataReader.h"
#include "itkQuadEdgeMeshScalarDataVTKPolyDataWriter.h"
#include "itkQuadEdgeMesh.h"

#include "SmoothingMeshScalarsCLP.h"
#include <BRAINSCommonLib.h>

int main( int argc, char *argv[] )
{
  PARSE_ARGS;
  BRAINSRegisterAlternateIO();

  if( inputSurfaceFile == "" )
    {
    std::cerr << "No input file specified" << std::endl;
    return 1;
    }
  if( outputSurfaceFile == "" )
    {
    std::cerr << "No output file specified" << std::endl;
    return 1;
    }

  std::cout << "---------------------------------------------------" << std::endl;
  std::cout << "Input Surface: " << inputSurfaceFile << std::endl;
  std::cout << "Output Surface: " << outputSurfaceFile << std::endl;
  std::cout << "Smooth the scalar values with lambda: " << lambda << std::endl;
  std::cout << "And the number of iterations: " << iterations << std::endl;
  std::cout << "---------------------------------------------------" << std::endl;

  const unsigned int Dimension = 3;
  typedef float MeshPixelType;

  typedef itk::QuadEdgeMesh<MeshPixelType, Dimension> InputMeshType;
  typedef itk::QuadEdgeMesh<MeshPixelType, Dimension> OutputMeshType;

  typedef itk::QuadEdgeMeshVTKPolyDataReader<InputMeshType> ReaderType;

  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( inputSurfaceFile.c_str() );

  typedef itk::QuadEdgeMeshScalarPixelValuesSmoothingFilter<
      InputMeshType, OutputMeshType>                       FilterType;

  FilterType::Pointer filter = FilterType::New();

  filter->SetLambda( lambda );
  filter->SetMaximumNumberOfIterations( iterations );

  typedef itk::QuadEdgeMeshScalarDataVTKPolyDataWriter<OutputMeshType> WriterType;
  WriterType::Pointer writer = WriterType::New();

  filter->SetInput( reader->GetOutput() );
  writer->SetInput( filter->GetOutput() );

  filter->Update();

  writer->SetFileName( outputSurfaceFile.c_str() );

  try
    {
    writer->Update();
    }
  catch( itk::ExceptionObject & excp )
    {
    std::cerr << excp << std::endl;
    return 1;
    }

  return 0;
}
