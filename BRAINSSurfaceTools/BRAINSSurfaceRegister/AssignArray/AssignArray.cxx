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
#include "vtkDataArray.h"
#include "vtkPointData.h"
#include "AssignArrayCLP.h"
#include "vtkSmartPointer.h"
#include <BRAINSCommonLib.h>

int main( int argc, char * argv[] )
{
  PARSE_ARGS;
  BRAINSRegisterAlternateIO();
  if( sourceSurfaceFile == "" )
    {
    std::cerr << "No source file specified" << std::endl;
    return 1;
    }
  if( targetSurfaceFile == "" )
    {
    std::cerr << "No target file specified" << std::endl;
    return 1;
    }
  if( outputSurfaceFile == "" )
    {
    std::cerr << "No output file specified" << std::endl;
    return 1;
    }

  if( getArrayName == "scalars" )
    {
    std::cout << "----------------------------------------------------" << std::endl;
    std::cout << "Assign scalars from: " << std::endl;
    std::cout << sourceSurfaceFile << std::endl;
    std::cout << "to: " << std::endl;
    std::cout << targetSurfaceFile << std::endl;
    std::cout << "With name of " << setArrayName << std::endl;
    std::cout << "Output File: " << std::endl;
    std::cout << outputSurfaceFile << std::endl;
    std::cout << "----------------------------------------------------" << std::endl;
    }
  else
    {
    std::cout << "----------------------------------------------------" << std::endl;
    std::cout << "Assign array: " << getArrayName << " from: " << std::endl;
    std::cout << sourceSurfaceFile << std::endl;
    std::cout << "to: " << std::endl;
    std::cout << targetSurfaceFile << std::endl;
    std::cout << "with name of " << setArrayName << std::endl;
    std::cout << "Output File: " << std::endl;
    std::cout << outputSurfaceFile << std::endl;
    std::cout << "----------------------------------------------------" << std::endl;
    }
  vtkSmartPointer<vtkPolyDataReader> s_reader = vtkSmartPointer<vtkPolyDataReader>::New();
  s_reader->SetFileName(sourceSurfaceFile.c_str() );
  s_reader->Update();

  vtkPolyData *source = s_reader->GetOutput();

  vtkSmartPointer<vtkPolyDataReader> t_reader = vtkSmartPointer<vtkPolyDataReader>::New();
  t_reader->SetFileName(targetSurfaceFile.c_str() );
  t_reader->Update();

  vtkPolyData *target = t_reader->GetOutput();

  vtkDataArray *label;

  if( getArrayName != "scalars" )
    {
    label = source->GetPointData()->GetArray(getArrayName.c_str() );
    }
  else
    {
    label = source->GetPointData()->GetScalars();
    }

  if( label == NULL )
    {
    std::cout << "there is no label array in the source." << std::endl;
    return 1;
    }

  // set new name to the label
  label->SetName(setArrayName.c_str() );

  if( target->GetPointData()->GetScalars() == NULL )
    {
    target->GetPointData()->SetScalars(label);
    }
  else
    {
    target->GetPointData()->AddArray(label);
    }

  vtkSmartPointer<vtkPolyDataWriter> writer = vtkSmartPointer<vtkPolyDataWriter>::New();
  writer->SetInput(target);
  writer->SetFileName(outputSurfaceFile.c_str() );
  writer->Update();

  return 0;
}
