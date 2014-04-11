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

#include "itkQuadEdgeMesh.h"
#include "itkIdentityTransform.h"

#include "itkLinearInterpolateMeshFunction.h"
#include "itkNearestNeighborInterpolateMeshFunction.h"
#include "itkResampleQuadEdgeMeshFilter.h"
#include "itkRescaleScalarsQuadEdgeMeshFilter.h"

#include "itkQuadEdgeMeshVTKPolyDataReader.h"
#include "itkQuadEdgeMeshScalarDataVTKPolyDataWriter.h"

#include "ResampleQuadEdgeMeshCLP.h"
#include <BRAINSCommonLib.h>

int main( int argc, char * argv [] )
{
  PARSE_ARGS;
  BRAINSRegisterAlternateIO();

  std::cout << "-----------------------------------------------" << std::endl;
  std::cout << "Input Mesh: " << std::endl;
  std::cout << inputMeshFile << std::endl;
  std::cout << "Reference Mesh: " << std::endl;
  std::cout << refMeshFile << std::endl;
  std::cout << "Output Mesh: " << std::endl;
  std::cout << outputMeshFile << std::endl;
  std::cout << "Interpolation Type: " << interpolateType << std::endl;
  std::cout << "-----------------------------------------------" << std::endl;

  typedef double PixelType;
  const unsigned int Dimension = 3;

  // Declaration of the type of Mesh
  typedef itk::QuadEdgeMesh<PixelType, Dimension> MeshType;

  // Here read a mesh from a file
  typedef itk::QuadEdgeMeshVTKPolyDataReader<MeshType> ReaderType;
  ReaderType::Pointer inputReader = ReaderType::New();
  inputReader->SetFileName( inputMeshFile.c_str() );
  inputReader->Update();

  ReaderType::Pointer refReader = ReaderType::New();
  refReader->SetFileName( refMeshFile.c_str() );
  refReader->Update();

  typedef itk::IdentityTransform<double, Dimension> TransformType;
  TransformType::Pointer transform = TransformType::New();
  transform->SetIdentity();

  typedef itk::LinearInterpolateMeshFunction<MeshType> LinearInterpolatorType;
  LinearInterpolatorType::Pointer interpolator_l = LinearInterpolatorType::New();

  typedef itk::NearestNeighborInterpolateMeshFunction<MeshType> NearestInterpolatorType;
  NearestInterpolatorType::Pointer interpolator_n = NearestInterpolatorType::New();

  typedef itk::ResampleQuadEdgeMeshFilter<MeshType, MeshType> MeshFilterType;
  MeshFilterType::Pointer filter = MeshFilterType::New();
  filter->SetReferenceMesh(refReader->GetOutput() );
  filter->SetInput(inputReader->GetOutput() );
  filter->SetTransform( transform );

  // set the interpolation type
  if( interpolateType == "Nearest" )
    {
    filter->SetInterpolator( interpolator_n );
    }
  else if( interpolateType == "Linear" )
    {
    filter->SetInterpolator( interpolator_l );
    }

  filter->Update();

  typedef itk::QuadEdgeMeshScalarDataVTKPolyDataWriter<MeshType> WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetInput( filter->GetOutput() );
  writer->SetFileName( outputMeshFile.c_str() );
  writer->Update();

  return 0;
}
