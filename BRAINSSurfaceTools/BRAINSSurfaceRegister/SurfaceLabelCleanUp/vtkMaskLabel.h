/*=========================================================================

  Program:   Visualization Toolkit
  Module:    $RCSfile: vtkMaskLabel.h,v $

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// .NAME vtkMaskLabel - keep subset of cells whose point has the specific
//                      label with it.
// .SECTION Description
// vtkMaskLabel is a filter that keeps the subset of cells if all of
// its corners have the specific label. The output will have
// a partition of the input cells but the whole members of points with
// point data.
//
// It requests label value
// has been saved as SCALARS in point data.
//
// Modified by Wen, 01/06/2012
// To add a flag to keep the output to have
// Only Labeled Points or Cells anly point of whose is Labeled

// .SECTION See Also
// vtkMaskPoints

#ifndef __vtkMaskLabel_h
#define __vtkMaskLabel_h

#include "vtkPolyDataAlgorithm.h"

class VTK_GRAPHICS_EXPORT vtkMaskLabel : public vtkPolyDataAlgorithm
{
public:
  static vtkMaskLabel * New();

  vtkTypeRevisionMacro(vtkMaskLabel, vtkPolyDataAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent);

  // Description:
  // Specify the label value
  vtkSetMacro(Label, int);
  vtkGetMacro(Label, int);

  // Turn on MaskPoint to keep only the "Labeled" points
  vtkSetMacro(LabelOnly, int);
  vtkGetMacro(LabelOnly, int);
  vtkBooleanMacro(LabelOnly, int);
protected:
  vtkMaskLabel();
  ~vtkMaskLabel()
  {
  };

  int RequestData(vtkInformation *, vtkInformationVector * *, vtkInformationVector *);

  int Label; // the label value that decides which cell is going to be kept.

  int LabelOnly; // boolean turns on/off to keep labeled points only
private:
  vtkMaskLabel(const vtkMaskLabel &);   // Not implemented.
  void operator=(const vtkMaskLabel &); // Not implemented.
};

#endif
