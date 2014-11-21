/*=========================================================================

Program:   Visualization Toolkit
Module:    vtkImageGenus0MarchingCubes.h

Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
All rights reserved.
See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

This software is distributed WITHOUT ANY WARRANTY; without even
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// .NAME vtkImageGenus0MarchingCubes - generate isosurface(s) from volume/images
// .SECTION Description
// vtkImageGenus0MarchingCubes is a filter that takes as input images (e.g., 3D
// image region) and generates on output one or more isosurfaces.
// One or more contour values must be specified to generate the isosurfaces.
// Alternatively, you can specify a min/max scalar range and the number of
// contours to generate a series of evenly spaced contour values.
// This filter can stream, so that the entire volume need not be loaded at
// once.  Streaming is controlled using the instance variable
// InputMemoryLimit, which has units KBytes.

// .SECTION Caveats
// This filter is specialized to volumes. If you are interested in
// contouring other types of data, use the general vtkContourFilter. If you
// want to contour an image (i.e., a volume slice), use vtkMarchingSquares.
// .SECTION See Also
// vtkContourFilter vtkSliceCubes vtkMarchingSquares vtkSynchronizedTemplates3D

// This is a modification fo the vtkImageMarchingCubes algorithm

#ifndef __vtkImageGenus0MarchingCubes_h
#define __vtkImageGenus0MarchingCubes_h

#include "vtkFiltersGeneralModule.h" // For export macro
#include "vtkPolyDataAlgorithm.h"
// -- #include <vtkTransform.h>

#include "vtkContourValues.h" // Needed for direct access to ContourValues

class vtkCellArray;
class vtkFloatArray;
class vtkImageData;
class vtkPoints;

#define VTK_BRAINSTOOLS_GRAPHICS_EXPORT
class VTK_BRAINSTOOLS_GRAPHICS_EXPORT vtkImageGenus0MarchingCubes : public vtkPolyDataAlgorithm
{
public:
  static vtkImageGenus0MarchingCubes * New();

  vtkTypeMacro(vtkImageGenus0MarchingCubes, vtkPolyDataAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent);

  // Description:
  // Methods to set contour values
  // why don't these do anything?
  void SetValue(int, double )
  {
  }

  double GetValue(int )
  {
    return 0.0;
  }

  int GetNumberOfConnectedComponents()
  {
    return iConnectedComponents;
  }

  vtkCellArray *Triangles;
  vtkPoints *   Points;

  vtkSetMacro(BiggestComponent, int );
  vtkGetMacro(BiggestComponent, int );
  vtkBooleanMacro(BiggestComponent, int );

  vtkSetMacro(ConnectedComponent, int );
  vtkGetMacro(ConnectedComponent, int );
  vtkBooleanMacro(ConnectedComponent, int );

  vtkSetMacro( CutLoops, int );
  vtkGetMacro( CutLoops, int );
  vtkBooleanMacro( CutLoops, int );

  void SetAltValue( int value )
  {
    altValue = value;
  }

  double GetAltValue()
  {
    return altValue;
  }

  vtkSetMacro( ComputeSurface, int );
  vtkGetMacro( ComputeSurface, int );
  vtkBooleanMacro( ComputeSurface, int );

  vtkSetMacro( Verbose, int );
  vtkGetMacro( Verbose, int );
  vtkBooleanMacro( Verbose, int );

  void Use18Connectivity()
  {
    iConnectivity = 18;
  }

  void Use6Connectivity()
  {
    iConnectivity = 6;
  }

  vtkImageData * GetCorrectedImageData()
  {
    return pCorrectedImageData;
  }

protected:
  vtkImageGenus0MarchingCubes();
  ~vtkImageGenus0MarchingCubes();

  void Execute();

  virtual int FillInputPortInformation(int port, vtkInformation *info);

  int BiggestComponent;
  int ConnectedComponent;
  int CutLoops;
  int Verbose;
  int ComputeSurface;

  int altValue;

  int iConnectivity;
  int iConnectedComponents;

  vtkImageData *pCorrectedImageData;
private:
  vtkImageGenus0MarchingCubes(const vtkImageGenus0MarchingCubes &); // Not
                                                                    // implemented.
  void operator=(const vtkImageGenus0MarchingCubes &);              // Not
                                                                    // implemented.

};

#endif
