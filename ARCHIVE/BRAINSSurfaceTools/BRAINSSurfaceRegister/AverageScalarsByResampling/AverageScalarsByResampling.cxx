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
#include "itkIdentityTransform.h"
#include "itkLinearInterpolateMeshFunction.h"
#include "itkResampleQuadEdgeMeshFilter.h"
#include "itkQuadEdgeMeshAddScalarsFilter.h"
#include "itkAssignScalarValuesQuadEdgeMeshFilter.h"

#include "AverageScalarsByResamplingCLP.h"
#include <BRAINSCommonLib.h>

int
main(int argc, char * argv[])
{
  PARSE_ARGS;
  BRAINSRegisterAlternateIO();
  // check the deformed template files
  if (int(deformedTemplateMeshList.size()) != numSubs)
  {
    std::cerr << "number of deformed meshes do not agree with number of subjects" << std::endl;
    return EXIT_FAILURE;
  }
  // check the sphere files with scalar values
  if (int(sphereWithScalarsList.size()) != numSubs)
  {
    std::cerr << "number of spheres with scalars do not agree with number of subjects" << std::endl;
    return EXIT_FAILURE;
  }
  // check template surface
  if (templateSurfaceFile == "")
  {
    std::cerr << "the template surface file should be specified" << std::endl;
    return EXIT_FAILURE;
  }
  // check the template sphere with scalars
  if (templateSphereFile == "")
  {
    std::cerr << "the template sphere with scalar values should be specified" << std::endl;
    return EXIT_FAILURE;
  }

  std::cout << "---------------------------------------------------" << std::endl;
  std::cout << "Add scalar values from: " << templateSphereFile << std::endl;
  std::cout << "Template sphere: " << templateSphereFile << std::endl;
  std::cout << "With scalar values from " << numSubs << " subjects calculated as: " << std::endl;
  for (int i = 0; i < numSubs; i++)
  {
    std::cout << "Use deformed template sphere: " << deformedTemplateMeshList[i] << std::endl;
    std::cout << "To resample scalar values of: " << sphereWithScalarsList[i] << std::endl;
  }
  std::cout << "Map the average scalars onto: " << templateSurfaceFile << std::endl;
  std::cout << "Output surface: " << templateSurfaceWithAverageScalars << std::endl;
  std::cout << "---------------------------------------------------" << std::endl;

  using MeshPixelType = double;
  constexpr unsigned int Dimension = 3;

  using MeshType = itk::QuadEdgeMesh<MeshPixelType, Dimension>;

  using InputMeshReaderType = itk::QuadEdgeMeshVTKPolyDataReader<MeshType>;

  // read input fixed mesh
  InputMeshReaderType::Pointer fixedMeshReader = InputMeshReaderType::New();
  fixedMeshReader->SetFileName(templateSurfaceFile.c_str());
  fixedMeshReader->Update();

  // set up deformed fixed mesh reader
  InputMeshReaderType::Pointer deformedFixedMeshReader = InputMeshReaderType::New();

  // set up reference mesh reader
  InputMeshReaderType::Pointer referenceMeshReader = InputMeshReaderType::New();

  // set up Interpolators Filter
  using TransformType = itk::IdentityTransform<double>;

  TransformType::Pointer transform = TransformType::New();

  using LinearInterpolatorType = itk::LinearInterpolateMeshFunction<MeshType>;

  LinearInterpolatorType::Pointer interpolator = LinearInterpolatorType::New();

  // set up Reample filter
  using ResamplingFilterType = itk::ResampleQuadEdgeMeshFilter<MeshType, MeshType>;

  ResamplingFilterType::Pointer resampleFilter = ResamplingFilterType::New();

  resampleFilter->SetTransform(transform);
  resampleFilter->SetInterpolator(interpolator);

  // set up Add Scalar Filter
  using AddScalarsFilterType = itk::QuadEdgeMeshAddScalarsFilter<MeshType, MeshType, MeshType>;

  AddScalarsFilterType::Pointer addFilter = AddScalarsFilterType::New();

  // set up Assign Scalar Filter
  using AssignFilterType = itk::AssignScalarValuesQuadEdgeMeshFilter<MeshType, MeshType, MeshType>;

  AssignFilterType::Pointer assignFilter = AssignFilterType::New();

  // read one deformed fixed mesh and one reference mesh once at a time
  // perform resampling and adding in one loop
  MeshType::Pointer deformedFixedMesh = MeshType::New();
  MeshType::Pointer referenceMesh = MeshType::New();
  MeshType::Pointer resampledMesh = MeshType::New();

  MeshType::Pointer newMesh = MeshType::New();

  InputMeshReaderType::Pointer inputWithScalarsReader = InputMeshReaderType::New();
  inputWithScalarsReader->SetFileName(templateSphereFile.c_str());
  inputWithScalarsReader->Update();

  MeshType::Pointer outputMesh = inputWithScalarsReader->GetOutput();
  outputMesh->DisconnectPipeline();
  for (int i = 0; i < numSubs; i++)
  {
    // read deformed fixed mesh
    deformedFixedMeshReader->SetFileName(deformedTemplateMeshList[i]);
    deformedFixedMeshReader->Update();

    deformedFixedMesh = deformedFixedMeshReader->GetOutput();
    deformedFixedMesh->DisconnectPipeline();

    // read reference mesh (sphere with scalar values)
    referenceMeshReader->SetFileName(sphereWithScalarsList[i]);
    referenceMeshReader->Update();

    referenceMesh = referenceMeshReader->GetOutput();
    referenceMesh->DisconnectPipeline();

    // resample reference mesh by deformed fixed mesh
    // resample the input by the reference
    resampleFilter->SetInput(referenceMesh);
    resampleFilter->SetReferenceMesh(deformedFixedMesh);
    resampleFilter->Update();

    resampledMesh = resampleFilter->GetOutput();
    resampledMesh->DisconnectPipeline();

    // assign the resampled scalars to input fixed mesh
    assignFilter->SetInputMesh(fixedMeshReader->GetOutput());
    assignFilter->SetSourceMesh(resampledMesh);
    assignFilter->Update();

    newMesh = assignFilter->GetOutput();
    newMesh->DisconnectPipeline();

    // add point.Value of newDF to point.Value of output
    addFilter->SetInput1(newMesh);
    addFilter->SetInput2(outputMesh);
    addFilter->Update();

    outputMesh = addFilter->GetOutput();

    outputMesh->DisconnectPipeline();
  }

  // divide the sum of scalars with numSubs
  // iterate through pointData of output
  using ScalarContainer = MeshType::PointDataContainer;

  ScalarContainer * scalars = outputMesh->GetPointData();

  using scalarIterator = ScalarContainer::Iterator;

  scalarIterator scalarItr = scalars->Begin();
  scalarIterator scalarEnd = scalars->End();

  while (scalarItr != scalarEnd)
  {
    if (numSubs != 0)
    {
      scalarItr.Value() = scalarItr.Value() / double(numSubs + 1);
    }
    ++scalarItr;
  }

  // write out deformation field
  using WriterType = itk::QuadEdgeMeshScalarDataVTKPolyDataWriter<MeshType>;
  WriterType::Pointer writer = WriterType::New();

  writer->SetInput(outputMesh);
  writer->SetFileName(templateSurfaceWithAverageScalars.c_str());
  writer->Update();

  return EXIT_SUCCESS;
}
