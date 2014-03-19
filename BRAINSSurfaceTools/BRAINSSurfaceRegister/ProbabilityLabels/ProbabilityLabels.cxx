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

#include "vtkPolyDataReader.h"
#include "vtkPolyDataWriter.h"
#include "vtkPolyData.h"
#include "vtkSmartPointer.h"
#include "vtkPointData.h"
#include "vtkDenseArray.h"
#include "vtkIntArray.h"
#include "vtkFloatArray.h"

#include "ProbabilityLabelsCLP.h"
#include <BRAINSCommonLib.h>

int main( int argc, char * argv[] )
{
  PARSE_ARGS;
  BRAINSRegisterAlternateIO();
  std::cout << "------------------------------------------" << std::endl;
  std::cout << "Input Meshes are: " << std::endl;
  for( unsigned int i = 0; i < inputMeshList.size(); i++ )
    {
    std::cout << inputMeshList[i] << std::endl;
    }
  std::cout << outputMeshFile << std::endl;
  if( mostLikely )
    {
    std::cout << "Calculate the most likely labels." << std::endl;
    }
  if( probability )
    {
    std::cout << "Calculate the label probabilities." << std::endl;
    }
  std::cout << "------------------------------------------" << std::endl;

  // read the first surface to get the idea about
  // label array, number of points, number of labels
  vtkSmartPointer<vtkPolyDataReader> polyDataReader = vtkSmartPointer<vtkPolyDataReader>::New();
  polyDataReader->SetFileName(inputMeshList[0].c_str() );
  polyDataReader->Update();

  // save it to templateMesh
  vtkSmartPointer<vtkPolyData> templateMesh = vtkSmartPointer<vtkPolyData>::New();
  templateMesh->DeepCopy(polyDataReader->GetOutput() );

  int    numOfPoints = templateMesh->GetNumberOfPoints();
  double labelRange[2];
  templateMesh->GetPointData()->GetScalars()->GetRange(labelRange);
  int startLabel = int(labelRange[0]);
  int endLabel = int(labelRange[1]);
  int numOfLabels = endLabel - startLabel + 1;

  // number of surfaces
  int numOfMeshes = inputMeshList.size();
  std::cout << "number of meshes: " << numOfMeshes << std::endl;

  // to keep all of the labels
  // using a 2D array of the size: numberOfSurfaces x numberOfPoints
  vtkSmartPointer<vtkDenseArray<float> > vertexLabels = vtkSmartPointer<vtkDenseArray<float> >::New();
  vertexLabels->Resize(numOfMeshes, numOfPoints);
  // read surfaces one by one
  // and keep label arrays in vertexLabels (size numberOfSurfaces x numberOfPoints)
  for( int i = 0; i < numOfMeshes; i++ )
    {
    // read the meshes
    // include the first one (template)
    polyDataReader->SetFileName(inputMeshList[i].c_str() );
    polyDataReader->Update();

    int nPoints = polyDataReader->GetOutput()->GetNumberOfPoints();

    if( nPoints != numOfPoints )
      {
      std::cerr << "Error: Invalid number of points in mesh: " << argv[i + 1];
      std::cerr << std::endl;
      std::cerr << "  Expected " << numOfPoints << "but found ";
      std::cerr << nPoints << std::endl;
      return 1;
      }
    // keep labels
    for( int j = 0; j < numOfPoints; j++ )
      {
      vertexLabels->SetValue(i, j, polyDataReader->GetOutput()->GetPointData()->GetScalars()->GetTuple1(j) );
      }
    }

  // if mostLikely is ON

  if( mostLikely )
    {
    // look at each column of vertexLabels
    // calculate frequency of each label between [startLabel,endLabel]
    // to save frequencies of labels at each point
    int *frequency = new int[numOfLabels];

    // create a new array to keep the most likely labels
    vtkSmartPointer<vtkIntArray> mostLikelyArray = vtkSmartPointer<vtkIntArray>::New();
    mostLikelyArray->SetNumberOfComponents(1);
    mostLikelyArray->SetNumberOfTuples(numOfPoints);
    mostLikelyArray->SetName("MostLikelyLabel");

    // create a new array to keep the number of different labels
    // happen at each point
    vtkSmartPointer<vtkIntArray> labelNumArray = vtkSmartPointer<vtkIntArray>::New();
    labelNumArray->SetNumberOfComponents(1);
    labelNumArray->SetNumberOfTuples(numOfPoints);
    labelNumArray->SetName("NumOfLabels");
    // go throught all meshes for one point
    for( int j = 0; j < numOfPoints; j++ )
      {
      // clean up frequency array
      for( int ii = 0; ii < numOfLabels; ii++ )
        {
        frequency[ii] = 0;
        }
      // count the frequency across meshes
      for( int i = 0; i < numOfMeshes; i++ )
        {
        int label_ij = vertexLabels->GetValue(i, j);
        frequency[label_ij - startLabel] += 1;
        }
      // search for the most likely label for point j
      int label_m = -1;
      int max_count = -1;
      int label_count = 0;
      for( int ii = 0; ii < numOfLabels; ii++ )
        {
        if( frequency[ii] != 0 )
          {
          label_count += 1;
          }
        if( frequency[ii] > max_count )
          {
          label_m = ii + startLabel;
          max_count = frequency[ii];
          }
        }
      // keep the most likely label to mostLikelyArray
      mostLikelyArray->SetTuple1(j, label_m);
      labelNumArray->SetTuple1(j, label_count);
      }

    // set the mostLikelyArray to be the scalars
    // overwrite what the template had
    templateMesh->GetPointData()->SetScalars(mostLikelyArray);
    templateMesh->GetPointData()->AddArray(labelNumArray);
    }
  // if probability is ON
  //               numOfSurfaces labeled as label_i
  // prabability = ---------------------------------
  //                   num of input surfaces
  for( int i = 0; i < numOfLabels; i++ )
    {
    int label_i = i + startLabel;
    // create a new array for label_i
    vtkSmartPointer<vtkFloatArray> labelPercentageArray = vtkSmartPointer<vtkFloatArray>::New();
    labelPercentageArray->SetNumberOfComponents(1);
    labelPercentageArray->SetNumberOfTuples(numOfPoints);
    // for each point, look at the probability of label_i
    for( int j = 0; j < numOfPoints; j++ )
      {
      int count = 0;
      for( int jj = 0; jj < numOfMeshes; jj++ )
        {
        if( vertexLabels->GetValue(jj, j) == label_i )
          {
          count += 1;
          }
        }
      // probability for point j, label i
      labelPercentageArray->SetTuple1(j, double(count) / double(numOfMeshes) );
      }

    // keep each label's probability using the label value
    char arrayName[20];
    sprintf(arrayName, "%d", i);
    labelPercentageArray->SetName(arrayName);

    templateMesh->GetPointData()->AddArray(labelPercentageArray);
    }
  // write outSurface
  vtkSmartPointer<vtkPolyDataWriter> writer = vtkSmartPointer<vtkPolyDataWriter>::New();
#if (VTK_MAJOR_VERSION < 6)
  writer->SetInput(templateMesh);
#else
  writer->SetInputData(templateMesh);
#endif
  writer->SetFileName(outputMeshFile.c_str() );
  writer->Update();

  return 0;
}
