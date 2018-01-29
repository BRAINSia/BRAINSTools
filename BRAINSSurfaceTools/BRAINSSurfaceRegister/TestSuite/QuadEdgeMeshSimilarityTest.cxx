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

This is a test for quad edge mesh similarity

 =========================================================================*/

#include "itkQuadEdgeMeshVTKPolyDataReader.h"
#include "itkQuadEdgeMeshSimilarityCalculator.h"
#include "itkQuadEdgeMesh.h"

int main( int argc, char * argv [] )
{
  if( argc < 5 )
    {
    std::cerr << "Missing arguments" << std::endl;
    std::cerr << "Usage: " << std::endl;
    std::cerr << argv[0] << std::endl;
    std::cerr << "inputMesh1 inputMesh2 LabelValue ExpectedDice";
    std::cerr << std::endl;
    return EXIT_FAILURE;
    }

  typedef int MeshPixelType;
  constexpr unsigned int Dimension = 3;

  typedef itk::QuadEdgeMesh<MeshPixelType, Dimension> InputMeshType1;
  typedef itk::QuadEdgeMesh<MeshPixelType, Dimension> InputMeshType2;

  typedef itk::QuadEdgeMeshVTKPolyDataReader<InputMeshType1> InputReaderType1;
  typedef itk::QuadEdgeMeshVTKPolyDataReader<InputMeshType2> InputReaderType2;

  InputReaderType1::Pointer inputMeshReader1 = InputReaderType1::New();
  inputMeshReader1->SetFileName( argv[1] );

  InputReaderType2::Pointer inputMeshReader2 = InputReaderType2::New();
  inputMeshReader2->SetFileName( argv[2] );

  try
    {
    inputMeshReader1->Update();
    inputMeshReader2->Update();
    }
  catch( itk::ExceptionObject & exp )
    {
    std::cerr << exp << std::endl;
    return EXIT_FAILURE;
    }

  typedef itk::QuadEdgeMeshSimilarityCalculator<InputMeshType1, InputMeshType2> SimilarityCalculatorType;

  SimilarityCalculatorType::Pointer similarityCalculator = SimilarityCalculatorType::New();

  similarityCalculator->SetInputMesh1(inputMeshReader1->GetOutput() );
  similarityCalculator->SetInputMesh2(inputMeshReader2->GetOutput() );

  similarityCalculator->SetLabelValue(atoi(argv[3]) );

  similarityCalculator->Compute();

  double diceTest = similarityCalculator->GetDice();
  double diceExp = atof(argv[4]);
  double tolerance = 1.0e-4;

  if( fabs(diceTest - diceExp) > tolerance )
    {
    std::cout << "Dice Test: " << diceTest << std::endl;
    std::cout << "Dice Expected: " << diceExp << std::endl;
    std::cout << "The resulted Dice is not expected." << std::endl;
    return EXIT_FAILURE;
    }

  return EXIT_SUCCESS;
}
