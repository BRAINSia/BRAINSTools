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
#include "vtkPointData.h"
#include "vtkSmartPointer.h"
#include "vtkDataArray.h"
#include "vtkVersionMacros.h"

#include "CombineLabelsCLP.h"
#include <BRAINSCommonLib.h>

int main( int argc, char * argv[] )
{
  PARSE_ARGS;
  BRAINSRegisterAlternateIO();

  // paramters:
  std::cout << "-----------------------------------------" << std::endl;
  std::cout << "Input Surface: " << inputSurfaceFile << std::endl;
  std::cout << "Output Surface: " << outputSurfaceFile << std::endl;
  std::cout << "Replace label: " << removeLabel << std::endl;
  std::cout << "With label: " << keepLabel << std::endl;
  std::cout << "-----------------------------------------" << std::endl;

  // read the input surface with prior label information
  vtkSmartPointer<vtkPolyDataReader> inputReader = vtkSmartPointer<vtkPolyDataReader>::New();
  inputReader->SetFileName(inputSurfaceFile.c_str() );
  inputReader->Update();

  vtkSmartPointer<vtkPolyData> inputSurface = inputReader->GetOutput();
  int                          nPoints = inputSurface->GetNumberOfPoints();

  // get label array from input surface
  vtkDataArray *labelArray = inputSurface->GetPointData()->GetArray("LabelValue");
  if( labelArray == NULL )
    {
    std::cerr << "There is no label array on the input surface. ";
    std::cerr << "Quit." << std::endl;
    return 1;
    }
  for( int i = 0; i < nPoints; i++ )
    {
    // if label_i == label2
    int label_i = labelArray->GetTuple1(i);
    if( label_i == removeLabel )
      {
      labelArray->SetTuple1(i, keepLabel);
      }
    }

  // write output Surface
  vtkSmartPointer<vtkPolyDataWriter> writer = vtkSmartPointer<vtkPolyDataWriter>::New();
#if (VTK_MAJOR_VERSION < 6)
  writer->SetInput(inputSurface);
#else
  writer->SetInputData(inputSurface);
#endif
  writer->SetFileName(outputSurfaceFile.c_str() );
  writer->Update();

  return 0;
}
