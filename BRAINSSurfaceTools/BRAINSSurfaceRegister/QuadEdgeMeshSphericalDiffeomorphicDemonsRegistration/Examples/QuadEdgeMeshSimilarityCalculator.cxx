/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkQuadEdgeMeshSimilarityCalculator.cxx,v $
  Language:  C++
  Date:      $Date: 2007-09-06 17:44:24 $
  Version:   $Revision: 1.3 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/

#ifdef _MSC_VER
#pragma warning ( disable : 4786 )
#endif

#include "itkQuadEdgeMeshVTKPolyDataReader.h"
#include "itkQuadEdgeMeshSimilarityCalculator.h"

#include "itkQuadEdgeMesh.h"

int main( int argc, char * argv [] )
{
  if( argc < 4 )
    {
    std::cerr << "Missing arguments" << std::endl;
    std::cerr << "Usage: " << std::endl;
    std::cerr << argv[0] << std::endl;
    std::cerr << "inputMesh1 inputMesh2 LabelValue";
    std::cerr << std::endl;
    return EXIT_FAILURE;
    }

  typedef int MeshPixelType;
  const unsigned int Dimension = 3;

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

  // print out the values of similarities
  std::cout << "the Dice measure of Label" << argv[3];
  std::cout << " is: " << similarityCalculator->GetDice() << std::endl;
  std::cout << "the Jaccard measure of Label" << argv[3];
  std::cout << " is: " << similarityCalculator->GetJaccard() << std::endl;

  return EXIT_SUCCESS;
}
