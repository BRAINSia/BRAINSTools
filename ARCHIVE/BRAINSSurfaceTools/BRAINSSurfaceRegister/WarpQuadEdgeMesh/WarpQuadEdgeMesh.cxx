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

#include "itkQuadEdgeMeshTraits.h"
#include "itkQuadEdgeMeshVTKPolyDataReader.h"
#include "itkQuadEdgeMeshScalarDataVTKPolyDataWriter.h"

#include "itkIdentityTransform.h"
#include "itkLinearInterpolateMeshFunction.h"
#include "itkNearestNeighborInterpolateMeshFunction.h"

#include "itkResampleQuadEdgeMeshFilter.h"
#include "itkAssignScalarValuesQuadEdgeMeshFilter.h"

#include "WarpQuadEdgeMeshCLP.h"
#include <BRAINSCommonLib.h>

int
main(int argc, char * argv[])
{
  PARSE_ARGS;
  BRAINSRegisterAlternateIO();

  std::cout << "-----------------------------------------------" << std::endl;
  std::cout << "Input Fixed Mesh: " << std::endl;
  std::cout << fixedMeshFile << std::endl;
  std::cout << "Input Moving Mesh: " << std::endl;
  std::cout << movingMeshFile << std::endl;
  std::cout << "Input Deformed Fixed Mesh: " << std::endl;
  std::cout << deformedMeshFile << std::endl;
  std::cout << "Output Mesh: " << std::endl;
  std::cout << outputMeshFile << std::endl;
  std::cout << "Interpolation Type: " << interpolateType << std::endl;
  std::cout << "-----------------------------------------------" << std::endl;

  using PixelType = float;
  constexpr unsigned int Dimension = 3;

  using MeshType = itk::QuadEdgeMesh<PixelType, Dimension>;

  using ReaderType = itk::QuadEdgeMeshVTKPolyDataReader<MeshType>;

  ReaderType::Pointer inputMeshReader = ReaderType::New();
  inputMeshReader->SetFileName(fixedMeshFile.c_str());
  inputMeshReader->Update();

  ReaderType::Pointer referenceMeshReader = ReaderType::New();
  referenceMeshReader->SetFileName(movingMeshFile.c_str());
  referenceMeshReader->Update();

  ReaderType::Pointer deformedMeshReader = ReaderType::New();
  deformedMeshReader->SetFileName(deformedMeshFile.c_str());
  deformedMeshReader->Update();

  using TransformType = itk::IdentityTransform<double>;

  TransformType::Pointer transform = TransformType::New();

  using LinearInterpolatorType = itk::LinearInterpolateMeshFunction<MeshType>;
  using NearestInterpolatorType = itk::NearestNeighborInterpolateMeshFunction<MeshType>;

  LinearInterpolatorType::Pointer interpolator_l = LinearInterpolatorType::New();

  NearestInterpolatorType::Pointer interpolator_n = NearestInterpolatorType::New();

  // get scalars from moving mesh (reference) for deformed mesh
  using ResamplingFilterType = itk::ResampleQuadEdgeMeshFilter<MeshType, MeshType>;

  ResamplingFilterType::Pointer resampler = ResamplingFilterType::New();

  resampler->SetTransform(transform);

  // set the interpolation type
  if (interpolateType == "Nearest")
  {
    resampler->SetInterpolator(interpolator_n);
  }
  else if (interpolateType == "Linear")
  {
    resampler->SetInterpolator(interpolator_l);
  }

  resampler->SetInput(referenceMeshReader->GetOutput());
  resampler->SetReferenceMesh(deformedMeshReader->GetOutput());

  resampler->Update();

  // assign scalars from deformed mesh to fixed mesh
  using AssignFilterType = itk::AssignScalarValuesQuadEdgeMeshFilter<MeshType, MeshType, MeshType>;

  AssignFilterType::Pointer assignFilter = AssignFilterType::New();

  assignFilter->SetInputMesh(inputMeshReader->GetOutput());
  assignFilter->SetSourceMesh(resampler->GetOutput());

  assignFilter->Update();

  // write the result
  using WriterType = itk::QuadEdgeMeshScalarDataVTKPolyDataWriter<MeshType>;
  WriterType::Pointer writer = WriterType::New();
  writer->SetInput(assignFilter->GetOutput());
  writer->SetFileName(outputMeshFile.c_str());
  writer->Update();

  return 0;
}
