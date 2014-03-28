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
#include "vtkDataArray.h"
#include "vtkStringArray.h"
#include "vtkFieldData.h"
#include "vtkVersionMacros.h"

#include <iostream> // I/O
#include <fstream>  // file I/O
#include <iomanip>  // format manipulation

#include "RearrangeSurfaceLabelsCLP.h"
#include <BRAINSCommonLib.h>

int main( int argc, char * argv[] )
{
  PARSE_ARGS;
  BRAINSRegisterAlternateIO();

  // read a surface with labels
  vtkSmartPointer<vtkPolyDataReader> reader = vtkSmartPointer<vtkPolyDataReader>::New();
  reader->SetFileName(inputSurfaceFile.c_str() );
  reader->Update();

  vtkSmartPointer<vtkPolyData> surface_in = reader->GetOutput();
  int                          nPoints = surface_in->GetNumberOfPoints();

  // check the labelArray
  vtkDataArray *labelArray = surface_in->GetPointData()->GetScalars();
  std::string   arrayName = labelArray->GetName();
  if( labelArray == NULL || arrayName != "LabelValue" )
    {
    std::cerr << "There is no label array exist on the surface. ";
    std::cerr << "Quit." << std::endl;
    return 1;
    }

  // calculate frequencies for labels between min and max
  double labelRange[2];
  labelArray->GetRange(labelRange);
  int  nLabels = int(labelRange[1] - labelRange[0] + 1);
  int *frequency = new int[nLabels];
  int *labelTransform = new int[nLabels];
  for( int i = 0; i < nLabels; i++ )
    {
    frequency[i] = 0;
    }
  for( int i = 0; i < nPoints; i++ )
    {
    int label_i = labelArray->GetTuple1(i);
    frequency[label_i] += 1;
    }

  // go through labels with non-zero frequencies
  int newLabel = 0;
  for( int i = 0; i < nLabels; i++ )
    {
    if( frequency[i] > 0 )
      {
      labelTransform[i] = newLabel;
      newLabel += 1;
      }
    }
  std::cout << "there are " << newLabel << " labels" << std::endl;
  // transform labelArray
  for( int i = 0; i < nPoints; i++ )
    {
    int label_i = labelArray->GetTuple1(i);
    labelArray->SetTuple1(i, labelTransform[label_i]);
    }

  vtkSmartPointer<vtkStringArray> nameArray = vtkSmartPointer<vtkStringArray>::New();
  if( labelNameFile != "" )
    {
    // attach label names from an external file
    // std::ifstream fin(labelNameFile);
    ifstream fin(labelNameFile.c_str() );
    // fin.open(labelNameFile.c_str(),std::ifstream::app);
    std::string labelName;

    nameArray->SetName("LabelName");
    nameArray->SetNumberOfValues(newLabel);
    for( int i = 0; i < newLabel; i++ )
      {
      // fin.getline(labelName, 100);
      getline(fin, labelName);
      if( labelName == "" )
        {
        std::cerr << "Number of label names does not match real number of labels." << std::endl;
        }
      else
        {
        nameArray->SetValue(i, labelName);
        std::cout << labelName << std::endl;
        }
      }
    }

  vtkSmartPointer<vtkFieldData> fd = vtkSmartPointer<vtkFieldData>::New();
  fd->AddArray(nameArray);

  surface_in->SetFieldData(fd);

  // write out the surface
  vtkSmartPointer<vtkPolyDataWriter> writer = vtkSmartPointer<vtkPolyDataWriter>::New();

#if (VTK_MAJOR_VERSION < 6)
  writer->SetInput(surface_in);
#else
  writer->SetInputData(surface_in);
#endif
  writer->SetFileName(outputSurfaceFile.c_str() );
  writer->Update();
  delete [] frequency;
  delete [] labelTransform;
  return 0;
}
