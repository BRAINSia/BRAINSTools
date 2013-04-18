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
#include "BRAINSvtkV6Compat.h"

#include "CombineLabelsCLP.h"

int main( int argc, char * argv[] )
{
  PARSE_ARGS;

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
  BRAINSvtkV6_SetInputData( writer, inputSurface);
  writer->SetFileName(outputSurfaceFile.c_str() );
  writer->Update();

  return 0;
}
