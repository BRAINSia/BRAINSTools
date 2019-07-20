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

#include "itkQuadEdgeMeshVTKPolyDataReader.h"
#include "itkQuadEdgeMeshScalarDataVTKPolyDataWriter.h"

#include "itkQuadEdgeMesh.h"
#include "itkQuadEdgeMeshClampScalarsFilter.h"

#include "QuadEdgeMeshClampScalarsCLP.h"
#include <BRAINSCommonLib.h>

int
main( int argc, char * argv[] )
{
  PARSE_ARGS;
  BRAINSRegisterAlternateIO();

  if ( inputSurfaceFile == "" )
  {
    std::cerr << "No input file specified" << std::endl;
    return 1;
  }
  if ( outputSurfaceFile == "" )
  {
    std::cerr << "No output file specified" << std::endl;
    return 1;
  }

  std::cout << "---------------------------------------------------" << std::endl;
  std::cout << "Input Surface: " << inputSurfaceFile << std::endl;
  std::cout << "Output Surface: " << outputSurfaceFile << std::endl;
  std::cout << "Clamp the scalar values into: " << std::endl;
  std::cout << "[ " << outputMin << " " << outputMax << " ]" << std::endl;
  std::cout << "---------------------------------------------------" << std::endl;

  using MeshPixelType = float;
  constexpr unsigned int Dimension = 3;

  using MeshType = itk::QuadEdgeMesh< MeshPixelType, Dimension >;

  using ReaderType = itk::QuadEdgeMeshVTKPolyDataReader< MeshType >;

  ReaderType::Pointer inputReader = ReaderType::New();
  inputReader->SetFileName( inputSurfaceFile.c_str() );
  inputReader->Update();

  using FilterType = itk::QuadEdgeMeshClampScalarsFilter< MeshType, MeshType >;

  FilterType::Pointer filter = FilterType::New();

  filter->SetInput( inputReader->GetOutput() );

  filter->ClampMinOn();
  filter->SetOutputMinimum( outputMin );

  filter->ClampMaxOn();
  filter->SetOutputMaximum( outputMax );

  filter->Update();

  using WriterType = itk::QuadEdgeMeshScalarDataVTKPolyDataWriter< MeshType >;
  WriterType::Pointer writer = WriterType::New();
  writer->SetInput( filter->GetOutput() );
  writer->SetFileName( outputSurfaceFile.c_str() );
  writer->Update();

  return 0;
}
