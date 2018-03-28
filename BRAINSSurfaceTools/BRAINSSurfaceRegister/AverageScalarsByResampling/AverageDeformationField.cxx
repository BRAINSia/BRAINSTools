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
// calculate the average deformation field of inputs
// inputs: deformation field as for the same fixed mesh
// output: the average deformation field on the same fixed mesh
// minimum parameters: NumberOfInputs DF1 output

#include "itkQuadEdgeMeshVTKPolyDataReader.h"
#include "itkQuadEdgeMeshVectorDataVTKPolyDataWriter.h"

#include "itkQuadEdgeMesh.h"
#include "itkQuadEdgeMeshTraits.h"
#include "itkQuadEdgeMeshAddDeformationFieldFilter.h"

int main( int argc, char * argv [] )
{
  if( argc < 4 )
    {
    std::cerr << "Missing arguments" << std::endl;
    std::cerr << "Usage: " << std::endl;
    std::cerr << argv[0] << std::endl;
    std::cerr << "NumberOfInputs DF1 DF2 ...";
    std::cerr << "outputDeformationField";
    std::cerr << std::endl;
    return EXIT_FAILURE;
    }

  unsigned int numInputs = atoi(argv[1]);

  typedef float MeshPixelType;
  constexpr unsigned int Dimension = 3;

  typedef itk::QuadEdgeMesh<MeshPixelType, Dimension> MeshType;

  typedef MeshType::PointType   PointType;
  typedef PointType::VectorType VectorType;

  typedef itk::QuadEdgeMeshTraits<VectorType, Dimension, bool, bool> VectorPointSetTraits;

  typedef itk::QuadEdgeMesh<VectorType, Dimension, VectorPointSetTraits> MeshWithVectorsType;

  typedef itk::QuadEdgeMeshVTKPolyDataReader<MeshWithVectorsType> DeformationFieldReaderType;

  DeformationFieldReaderType::Pointer deformationFieldReader = DeformationFieldReaderType::New();

  typedef itk::QuadEdgeMeshAddDeformationFieldFilter<
      MeshWithVectorsType, MeshWithVectorsType, MeshWithVectorsType> AddDeformationFieldFilterType;

  AddDeformationFieldFilterType::Pointer addDFFilter = AddDeformationFieldFilterType::New();

  // read one deformationField once in a time
  // add it to the previous deformationField
  MeshWithVectorsType::Pointer output = MeshWithVectorsType::New();
  MeshWithVectorsType::Pointer newDF = MeshWithVectorsType::New();
  for( unsigned int i = 0; i < numInputs; i++ )
    {
    deformationFieldReader->SetFileName( argv[i + 2] );
    deformationFieldReader->Update();

    std::cout << "read deformation field: ";
    std::cout << argv[i + 2] << std::endl;

    if( i == 0 )
      {
      output = deformationFieldReader->GetOutput();
      output->DisconnectPipeline();
      }
    else
      {
      newDF = deformationFieldReader->GetOutput();
      newDF->DisconnectPipeline();

      // add point.Value of newDF to point.Value of output
      addDFFilter->SetInput1(output);
      addDFFilter->SetInput2(newDF);
      addDFFilter->Update();

      output = addDFFilter->GetOutput();

      output->DisconnectPipeline();
      }
    }

  // divide the sum deformation field with numInputs
  // iterate through pointData of output
  typedef MeshWithVectorsType::PointDataContainer DisplacementVectorContainer;

  DisplacementVectorContainer * displacementVectors = output->GetPointData();

  typedef DisplacementVectorContainer::Iterator displacementVectorIterator;

  displacementVectorIterator displacementVectorItr = displacementVectors->Begin();
  displacementVectorIterator displacementVectorEnd = displacementVectors->End();

  while( displacementVectorItr != displacementVectorEnd )
    {
    if( numInputs != 0 )
      {
      displacementVectorItr.Value() = displacementVectorItr.Value() / double(numInputs);
      }
    ++displacementVectorItr;
    }

  // write out deformation field
  typedef itk::QuadEdgeMeshVectorDataVTKPolyDataWriter<MeshWithVectorsType> VectorMeshWriterType;
  VectorMeshWriterType::Pointer vectorMeshWriter = VectorMeshWriterType::New();

  vectorMeshWriter->SetInput( output );
  vectorMeshWriter->SetFileName(argv[numInputs + 2]);
  vectorMeshWriter->Update();

  return EXIT_SUCCESS;
}
