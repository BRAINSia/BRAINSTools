/*=========================================================================

  Program:   Visualization Toolkit
  Module:    $RCSfile: vtkStructuredGrid.h,v $

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
#ifndef __vtkTalairachGrid_h
#define __vtkTalairachGrid_h

#include "vtkPoints.h"
#include "vtkStructuredGrid.h"
#include <vector>

class VTK_FILTERING_EXPORT vtkTalairachGrid : public vtkStructuredGrid
{
public:
  static vtkTalairachGrid * New();

  vtkTypeRevisionMacro(vtkTalairachGrid, vtkDataSet);

  /* Description:
   * Clear out the memory associated with both the box and grid */
  void Initialize();

  /* Description:
   * Get the dimensions of the bounding box and the talairach grid */
  int * GetBoundingBoxDimensions();

  int * GetTalairachGridDimensions();

  /* Description:
   * Set the AC, PC, SLA, and IRP Points. These points define the talairach
   * grid. Call EstablishGrid() to create the underlying vtkStructuredGrid */
  void SetACPoint(double ac[3]);
  void SetPCPoint(double pc[3]);
  void SetIRPPoint(double irp[3]);
  void SetSLAPoint(double sla[3]);

  /* Description:
   * Get the AC, PC, SLA, and IRP Points. These points define the
   * talairach grid. */
  double * GetACPoint();

  double * GetPCPoint();

  double * GetIRPPoint();

  double * GetSLAPoint();

  /* Description:
   * Get the underlying vtkStructuredGrid representations for the Talairach
   * coordinate system */
  vtkStructuredGrid * GetTalairachGrid();

  vtkStructuredGrid * GetBoundingBoxGrid();

  /* Description:
   * Set the underlying vtkStructuredGrid representations for the Talairach
   * coordinate system */
  void SetTalairachGrid( vtkStructuredGrid * );

  void SetBoundingBoxGrid( vtkStructuredGrid * );

  /* Description:
   * Get the points list for the box and grid */
  vtkPoints * GetTalairachGridPoints();

  vtkPoints * GetBoundingBoxGridPoints();

  /* Description:
   * Insert the points for both the box and grid using the provided AC, PC,
   * SLA and IRP points */
  void EstablishTalairachGrid();

  void EstablishBoundingBoxGrid();

  /* Description:
   * Set the points and generate both the box and grid */
  void Update();

  /* Description:
   * Write out the box and grid to a file */
  void WriteTalairachGrid(std::string filename);

  void WriteBoundingBoxGrid(std::string filename);

  /* Description:
   * Write out the appropriate data when the talairachGrid object is added to
   * an IO stream */
  void PrintSelf(ostream & os, vtkIndent indent);

protected:

  vtkTalairachGrid();
  ~vtkTalairachGrid();
private:

  /* Convert between talairach and voxel points */
  std::vector<double> ConvertTalairachPointToPixelPoint(double *talPoint);

  std::vector<double> ConvertPixelPointToTalairachPoint(double *voxelPoint);

  /* Stores the 4 points used to define all other points in the box and grid */
  double ACPoint[3];
  double PCPoint[3];
  double SLAPoint[3];
  double IRPPoint[3];

  /* Stores the complete set of points that define the box and grid */
  vtkPoints *boundingBoxGridPoints;
  vtkPoints *talairachGridPoints;

  /* The two Grid representations */
  vtkStructuredGrid *boundingBoxGrid;
  vtkStructuredGrid *talairachGrid;

  /* The data extent */
  int Extent[6];

  vtkTalairachGrid(const vtkTalairachGrid &); /* Not implemented. */
  void operator=(const vtkTalairachGrid &);   /* Not implemented. */

};

#endif
