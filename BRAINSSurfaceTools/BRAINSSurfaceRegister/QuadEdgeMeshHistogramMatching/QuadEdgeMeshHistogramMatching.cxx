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
#include "itkHistogramMatchingQuadEdgeMeshFilter.h"
#include "QuadEdgeMeshHistogramMatchingCLP.h"
#include <BRAINSCommonLib.h>

int main( int argc, char * argv [] )
{
  PARSE_ARGS;
  BRAINSRegisterAlternateIO();

  if( inputSurfaceFile == "" )
    {
    std::cerr << "No input file specified" << std::endl;
    return 1;
    }
  if( refSurfaceFile == "" )
    {
    std::cerr << "No reference file specified" << std::endl;
    return 1;
    }
  if( outputSurfaceFile == "" )
    {
    std::cerr << "No output file specified" << std::endl;
    return 1;
    }

  std::cout << "---------------------------------------------------" << std::endl;
  std::cout << "Input Surface: " << inputSurfaceFile << std::endl;
  std::cout << "Reference Surface: " << refSurfaceFile << std::endl;
  std::cout << "Output Surface: " << outputSurfaceFile << std::endl;
  std::cout << "Adjust scalar values of the input surface ";
  std::cout << "according to the reference surface ";
  std::cout << "using a histogram matching technique with ";
  std::cout << "number of histogram levels of " << numberOfHistogramLevels;
  std::cout << " and number of match points of " << numberOfMatchPoints << std::endl;
  std::cout << "---------------------------------------------------" << std::endl;

  typedef float MeshPixelType;
  constexpr unsigned int Dimension = 3;

  typedef itk::QuadEdgeMesh<MeshPixelType, Dimension> MeshType;

  typedef itk::QuadEdgeMeshVTKPolyDataReader<MeshType> SourceReaderType;
  typedef itk::QuadEdgeMeshVTKPolyDataReader<MeshType> ReferenceReaderType;

  SourceReaderType::Pointer srcReader = SourceReaderType::New();
  srcReader->SetFileName( inputSurfaceFile.c_str() );

  ReferenceReaderType::Pointer refReader = ReferenceReaderType::New();
  refReader->SetFileName( refSurfaceFile.c_str() );

  srcReader->Update();
  refReader->Update();

  typedef itk::HistogramMatchingQuadEdgeMeshFilter<MeshType, MeshType> FilterType;

  FilterType::Pointer filter  = FilterType::New();

  filter->SetInput(srcReader->GetOutput() );
  filter->SetReferenceMesh(refReader->GetOutput() );
  filter->SetNumberOfHistogramLevels( numberOfHistogramLevels );
  filter->SetNumberOfMatchPoints( numberOfMatchPoints );
  filter->Update();

  typedef itk::QuadEdgeMeshScalarDataVTKPolyDataWriter<MeshType> WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetInput( filter->GetOutput() );
  writer->SetFileName(outputSurfaceFile.c_str() );
  writer->Update();

  return 0;
}
