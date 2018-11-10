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

  Program:   Visualization Toolkit
  Module:    $RCSfile: vtkMaskLabel.cxx,v $

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
#include "vtkMaskLabel.h"

#include "vtkCellArray.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkObjectFactory.h"
#include "vtkPointData.h"
#include "vtkPolyData.h"

#include "itkMacro.h" //Needed for nullptr

vtkStandardNewMacro(vtkMaskLabel);

vtkMaskLabel::vtkMaskLabel()
{
  this->Label = 0;
  this->LabelOnly = 0;
}

// cut the cells. don't cut the points
// also keep all of the point data with them.
int vtkMaskLabel::RequestData(
  vtkInformation *vtkNotUsed(request),
  vtkInformationVector * *inputVector,
  vtkInformationVector *outputVector)
{
  // get the info objects
  vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
  vtkInformation *outInfo = outputVector->GetInformationObject(0);

  // get the input and output
  vtkPolyData *input = vtkPolyData::SafeDownCast(
      inInfo->Get(vtkDataObject::DATA_OBJECT() ) );
  vtkPolyData *output = vtkPolyData::SafeDownCast(
      outInfo->Get(vtkDataObject::DATA_OBJECT() ) );

  vtkIdType     id, pid;
  vtkPointData *pd;
  vtkIdType     numCells;
  vtkIdType *   pts = nullptr;
  vtkIdType     npts = 0;
  bool          LabelYes;
  int           abortExecute = 0;

  // Check input / pass data through
  numCells = input->GetNumberOfCells();

  if( numCells < 1 )
    {
    vtkErrorMacro(<< "No PolyData to mask!");
    return 1;
    }

  output->Allocate(input, numCells);
  input->BuildCells();

  // Traverse topological lists and traverse
  vtkIdType tenth = numCells / 10 + 1;
  for( id = 0; id < numCells && !abortExecute; id++ )
    {
    if( !(id % tenth) )
      {
      this->UpdateProgress( (float)id / numCells);
      abortExecute = this->GetAbortExecute();
      }
    input->GetCellPoints(id, npts, pts);

    // use LabelOnly to decide how to insert cells to output
    if( this->LabelOnly )
      {
      LabelYes = true;
      for( pid = 0; pid < npts; pid++ )
        {
        if( int(input->GetPointData()->GetScalars()->GetTuple1(pts[pid]) )
            != this->Label )
          {
          LabelYes = false;
          }
        }
      if( LabelYes )
        {
        output->InsertNextCell(input->GetCellType(id), npts, pts);
        }
      }
    else
      {
      LabelYes = false;
      for( pid = 0; pid < npts; pid++ )
        {
        if( int(input->GetPointData()->GetScalars()->GetTuple1(pts[pid]) )
            == this->Label )
          {
          LabelYes = true;
          }
        }
      if( LabelYes )
        {
        output->InsertNextCell(input->GetCellType(id), npts, pts);
        }
      }
    }

  // Update ourselves and release memory
  output->SetPoints(input->GetPoints() );
  pd = input->GetPointData();
  output->GetPointData()->PassData(pd);

  output->Squeeze();

  return 1;
}

void vtkMaskLabel::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);

  os << indent << "Label: " << this->Label << "\n";
}
