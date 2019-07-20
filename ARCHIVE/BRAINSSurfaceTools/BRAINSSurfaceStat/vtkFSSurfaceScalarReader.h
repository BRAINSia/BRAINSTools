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
/*=auto=========================================================================

  Portions (c) Copyright 2005 Brigham and Women's Hospital (BWH) All Rights Reserved.

  See Doc/copyright/copyright.txt
  or http://www.slicer.org/copyright/copyright.txt for details.

  Program:   3D Slicer
  Module:    $RCSfile: vtkFSSurfaceScalarReader.h,v $
  Date:      $Date: 2006/05/26 19:40:15 $
  Version:   $Revision: 1.9 $

=========================================================================auto=*/
/// .NAME vtkFSSurfaceScalarReader - read a surface scalar file from
/// Freesurfer tools
/// .SECTION Description
/// Reads a surface scalar file from FreeSurfer and outputs a
/// vtkFloatArray. Use the SetFileName function to specify the file
/// name. The number of values in the array should be equal to the
/// number of vertices/points in the surface.

#ifndef __vtkFSSurfaceScalarReader_h
#define __vtkFSSurfaceScalarReader_h

#include "vtkDataReader.h"
#include "vtkPolyData.h"
#include "vtkFloatArray.h"

class vtkFloatArray;
class vtkFSSurfaceScalarReader : public vtkDataReader
{
public:
  static vtkFSSurfaceScalarReader *
  New();

  vtkTypeMacro( vtkFSSurfaceScalarReader, vtkDataReader );
  void
  PrintSelf( ostream & os, vtkIndent indent );

  vtkFloatArray *
  GetOutput()
  {
    return this->scalars;
  };
  void
  SetOutput( vtkFloatArray * output )
  {
    this->scalars = output;
  };

  ///
  /// Read the scalars from a file. Return 1 on success, 0 on failure
  int
  ReadFSScalars();

  /// file type magic numbers
  /// const int FS_NEW_SCALAR_MAGIC_NUMBER = 16777215;
  // BTX
  enum
  {
    FS_NEW_SCALAR_MAGIC_NUMBER = 16777215,
  };
  // ETX
protected:
  vtkFSSurfaceScalarReader();
  ~vtkFSSurfaceScalarReader();

  vtkFloatArray * scalars;

  int
  ReadInt3( FILE * iFile, int & oInt );

  int
  ReadInt2( FILE * iFile, int & oInt );

  int
  ReadFloat( FILE * iFile, float & oInt );

private:
  vtkFSSurfaceScalarReader( const vtkFSSurfaceScalarReader & ); /// Not implemented.
  void
  operator=( const vtkFSSurfaceScalarReader & ); /// Not implemented.
};

#endif
