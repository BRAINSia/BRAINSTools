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
#include "itkConvertVTKToQuadEdgeMeshFilter.h"
#include "itkQuadEdgeMeshVTKPolyDataReader.h"
#include "itkQuadEdgeMeshSimilarityCalculator.h"
#include "itkQuadEdgeMesh.h"

#include "vtkPolyDataReader.h"
#include "vtkPolyData.h"
#include "vtkSmartPointer.h"
#include "vtkPointData.h"

#include <iostream> // I/O
#include <fstream>  // file I/O
#include <iomanip>  // format manipulation

#include "QuadEdgeMeshSimilarityCLP.h"
#include <BRAINSCommonLib.h>

int main( int argc, char * argv [] )
{
  PARSE_ARGS;
  BRAINSRegisterAlternateIO();
  std::cout << "--------------------------------------------" << std::endl;
  std::cout << "Caculate the Dice and Jaccard indices between " << std::endl;
  std::cout << labelName1 << " of surface: " << inputSurfaceFile1 << std::endl;
  std::cout << "And " << labelName2 << " of surface: " << inputSurfaceFile2 << std::endl;
  std::cout << "Output .txt file: " << outputSimilarityFile << std::endl;
  if( average )
    {
    std::cout << "Also output the average indices across labels." << std::endl;
    }
  std::cout << "--------------------------------------------" << std::endl;

  // use vtk to get the range of labels for both inputs
  // and take care of label arrays if they are not saved as scalars
  vtkSmartPointer<vtkPolyDataReader> vtkReader1 = vtkSmartPointer<vtkPolyDataReader>::New();
  vtkReader1->SetFileName(inputSurfaceFile1.c_str() );
  vtkReader1->Update();

  vtkPolyData * inputSurface1 = vtkReader1->GetOutput();
  int           numOfPoints = inputSurface1->GetNumberOfPoints();
  vtkDataArray *labelArray1;
  if( labelName1 == "scalars" )
    {
    labelArray1 = inputSurface1->GetPointData()->GetScalars();
    }
  else
    {
    labelArray1 = inputSurface1->GetPointData()->GetArray(labelName1.c_str() );
    }
  if( labelArray1 == NULL )
    {
    std::cerr << "surface1 does't have label array with the name: " << labelName1 << std::endl;
    std::cerr << "Quit." << std::endl;
    return 1;
    }
  double labelRange1[2];
  labelArray1->GetRange(labelRange1);

  vtkSmartPointer<vtkPolyDataReader> vtkReader2 = vtkSmartPointer<vtkPolyDataReader>::New();
  vtkReader2->SetFileName(inputSurfaceFile2.c_str() );
  vtkReader2->Update();

  vtkPolyData *inputSurface2 = vtkReader2->GetOutput();
  if( inputSurface2->GetNumberOfPoints() != numOfPoints )
    {
    std::cerr << "Number of points on surface2 is different from surface1." << std::endl;
    std::cerr << "Quit." << std::endl;
    return 1;
    }
  vtkDataArray *labelArray2;
  if( labelName2 == "scalars" )
    {
    labelArray2 = inputSurface2->GetPointData()->GetScalars();
    }
  else
    {
    labelArray2 = inputSurface2->GetPointData()->GetArray(labelName2.c_str() );
    }
  if( labelArray2 == NULL )
    {
    std::cerr << "surface2 does't have label array with the name: " << labelName2 << std::endl;
    std::cerr << "Quit." << std::endl;
    return 1;
    }
  double labelRange2[2];
  labelArray2->GetRange(labelRange2);
  if( labelRange1[0] != labelRange2[0] || labelRange1[1] != labelRange2[1] )
    {
    std::cerr << "The label range on surface2 is different from surface1." << std::endl;
    std::cerr << "Quit." << std::endl;
    return 1;
    }

  // set labelArray as scalars
  inputSurface1->GetPointData()->SetActiveScalars(labelName1.c_str() );
  inputSurface2->GetPointData()->SetActiveScalars(labelName2.c_str() );

  // convert vtk polydata to itk QuadEdgeMesh
  // use integer as pixelType of QuadEdgeMesh

  typedef int MeshPixelType;
  const unsigned int Dimension = 3;

  typedef itk::QuadEdgeMesh<MeshPixelType, Dimension> MeshType;

  typedef itk::ConvertVTKToQuadEdgeMeshFilter<MeshType> convertFilterType;

  convertFilterType::Pointer convertor1 = convertFilterType::New();
  convertor1->SetPolyData(inputSurface1);
  convertor1->Update();

  convertFilterType::Pointer convertor2 = convertFilterType::New();
  convertor2->SetPolyData(inputSurface2);
  convertor2->Update();

  // set up QuadEdgeMesh similarity calculator
  typedef itk::QuadEdgeMeshSimilarityCalculator<MeshType, MeshType> SimilarityCalculatorType;

  SimilarityCalculatorType::Pointer similarityCalculator = SimilarityCalculatorType::New();

  similarityCalculator->SetInputMesh1( convertor1->GetOutput() );
  similarityCalculator->SetInputMesh2( convertor2->GetOutput() );

  double sum_dice = 0.0;
  double sum_jaccard = 0.0;

  int startLabel = int(labelRange1[0]);
  int endLabel = int(labelRange1[1]);
  int numOfLabels = endLabel - startLabel + 1;

  // go through labels
  std::ofstream similarity_out;
  similarity_out.open(outputSimilarityFile.c_str(), std::ofstream::app);
  for( int i = startLabel; i < numOfLabels; i++ )
    {
    similarityCalculator->SetLabelValue(i);
    similarityCalculator->Compute();

    sum_dice += similarityCalculator->GetDice();
    sum_jaccard += similarityCalculator->GetJaccard();

    // write out the similarity file
    similarity_out << i << ": " << similarityCalculator->GetDice();
    similarity_out << " " << similarityCalculator->GetJaccard() << std::endl;
    }

  if( average )
    {
    similarity_out << "Average Dice: " << sum_dice / double(numOfLabels) << std::endl;
    similarity_out << "Average Jaccard: " << sum_jaccard / double(numOfLabels) << std::endl;
    }
  similarity_out.close();

  return 0;
}
