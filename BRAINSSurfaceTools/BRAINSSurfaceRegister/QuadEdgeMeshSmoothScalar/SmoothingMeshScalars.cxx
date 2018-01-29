/*=========================================================================
 *
 *  Copyright SINAPSE: Scalable Informatics for Neuroscience, Processing and Software Engineering
 *            The University of Iowa
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         http://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *=========================================================================*/
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

  constexpr unsigned int Dimension = 3;
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
