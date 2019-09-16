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
  Module:    $RCSfile: vtkFSSurfaceReader.h,v $
  Date:      $Date: 2006/05/26 19:40:14 $
  Version:   $Revision: 1.9 $

=========================================================================auto=*/
/// .NAME vtkFSSurfaceReader - read a surface file from Freesurfer tools
/// .SECTION Description
/// Reads a surface file from FreeSurfer and output PolyData. Use the
/// SetFileName function to specify the file name.

#ifndef __vtkFSSurfaceReader_h
#define __vtkFSSurfaceReader_h

#include "vtkDataReader.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkPolyData.h"

/// Prints debugging info.
#define FS_DEBUG 0

/// This code can calc normals but it doesn't seem to do so very well
/// (bug?) and vtkPolyDataNormals does it anyway.
#define FS_CALC_NORMALS 0

class vtkInformation;
class vtkInformationVector;
class vtkPolyData;
class vtkFSSurfaceReader : public vtkDataReader
{
public:
  static vtkFSSurfaceReader *
  New();

  vtkTypeMacro(vtkFSSurfaceReader, vtkDataReader);
  void
  PrintSelf(ostream & os, vtkIndent indent);

  ///
  /// Get the output of this reader.
  vtkPolyData *
  GetOutput();

  vtkPolyData *
  GetOutput(int idx);

  void
  SetOutput(vtkPolyData * output);

  // BTX
  // Description:
  /// old previous versions constants
  enum
  {
    FS_QUAD_FILE_MAGIC_NUMBER = (-1 & 0x00ffffff),
    FS_NEW_QUAD_FILE_MAGIC_NUMBER = (-3 & 0x00ffffff),
    FS_TRIANGLE_FILE_MAGIC_NUMBER = (-2 & 0x00ffffff),
    FS_NUM_VERTS_IN_QUAD_FACE = 4,    /// dealing with quads
    FS_NUM_VERTS_IN_TRI_FACE = 3,     /// dealing with tris
    FS_MAX_NUM_FACES_PER_VERTEX = 10, /// kinda arbitrary
  };
  // ETX

  int
  RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector * outputVector);

protected:
  vtkFSSurfaceReader();
  ~vtkFSSurfaceReader();

  // void Execute();

  /// Update extent of PolyData is specified in pieces.
  /// Since all DataObjects should be able to set UpdateExent as pieces,
  /// just copy output->UpdateExtent  all Inputs.
#if (VTK_MAJOR_VERSION >= 5)
  virtual int
  FillOutputPortInformation(int, vtkInformation *);

#else
  void
  ComputeInputUpdateExtents(vtkDataObject * output);

#endif

  /// Used by streaming: The extent of the output being processed by
  /// the execute method. Set in the ComputeInputUpdateExtents method.
  int ExecutePiece;
  int ExecuteNumberOfPieces;
  int ExecuteGhostLevel;

private:
  vtkFSSurfaceReader(const vtkFSSurfaceReader &); /// Not implemented.
  void
  operator=(const vtkFSSurfaceReader &); /// Not implemented.
};

/// If we're going to try to do normals, we need to keep some
/// connectivity information for the verts and faces. Use these
/// structures.
#if FS_CALC_NORMALS
typedef struct
{
  int   numFaces;
  int   faces[FS_MAX_NUM_FACES_PER_VERTEX];
  int   indicesInFace[FS_MAX_NUM_FACES_PER_VERTEX];
  float x, y, z;
  float nx, ny, nz;
} Vertex;

typedef struct
{
  int vertices[FS_NUM_SIDES_IN_FACE];
} Face;
#endif

#endif
