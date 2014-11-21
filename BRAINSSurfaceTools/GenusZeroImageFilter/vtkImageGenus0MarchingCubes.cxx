/*=========================================================================

  Program:   Visualization Toolkit
  Module:    $RCSfile: vtkImageGenus0MarchingCubes.cxx,v $

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
#include "vtkImageGenus0MarchingCubes.h"

#include "vtkCellArray.h"
#include "vtkCommand.h"
#include "vtkFloatArray.h"
#include "vtkImageData.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
// #include "vtkMarchingCubesCases.h"
#include "vtkObjectFactory.h"
#include "vtkPointData.h"
#include "vtkPolyData.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include <math.h>
#include <map>

#include "genus0.h"

vtkCxxRevisionMacro(vtkImageGenus0MarchingCubes, "$Revision: 1.1 $");
vtkStandardNewMacro(vtkImageGenus0MarchingCubes);

// ----------------------------------------------------------------------------
template <class T>
void vtkGetUnsignedShortData(vtkImageGenus0MarchingCubes *,
                             vtkImageData *, T *ptr,
                             unsigned short *pus, int totlen)
{
  for ( int iI = 0; iI < totlen; iI++ )
    {
    if ( ptr[iI] > 0 )
      {
      pus[iI] = 1;
      }
    else
      {
      pus[iI] = 0;
      }
    }
}

// ----------------------------------------------------------------------------
template <class T>
void vtkSetUnsignedShortData(vtkImageGenus0MarchingCubes *,
                             vtkImageData *, T *ptr,
                             unsigned short *pus, int totlen)
{
  for ( int iI = 0; iI < totlen; iI++ )
    {
    ptr[iI] = pus[iI];
    }
}

// ----------------------------------------------------------------------------
// Description:
// Construct object with initial range (0,1) and single contour value
// of 0.0. ComputeNormal is on, ComputeGradients is off and ComputeScalars is
// on.
vtkImageGenus0MarchingCubes::vtkImageGenus0MarchingCubes()
{
  BiggestComponent = 1;
  ConnectedComponent = 0;
  CutLoops = 0;
  Verbose = 1;
  ComputeSurface = 0;

  altValue = 1;

  // IJKtoRAS = NULL;

  iConnectivity = 18;

  iConnectedComponents = 0;

  pCorrectedImageData = NULL;
}

vtkImageGenus0MarchingCubes::~vtkImageGenus0MarchingCubes()
{
  if ( pCorrectedImageData != NULL )
    {
    pCorrectedImageData->Delete();
    pCorrectedImageData = NULL;
    }
}

void vtkImageGenus0MarchingCubes::Execute()
{
  vtkImageData *inData = vtkImageData::SafeDownCast( this->GetInput( 0 ) );
  vtkPolyData  *outData = vtkPolyData::SafeDownCast( this->GetOutput( 0 ) );

  this->iConnectedComponents = 0;

  // Null input check
  if ( !inData )
    {
    std::cerr << "Error: Input data not set." << std::endl;
    return;
    }

  // vtkMatrix4x4* matIJKtoRAS = IJKtoRAS->GetMatrix();

  // need to get the input dimensions and set up the Haker Genus 0 code

  int dims[3];
  inData->GetDimensions( dims );  // get the dimensions

  std::cout << "Dimensions = " << dims[0] << " " << dims[1] << " " << dims[2] << std::endl;

  double dSpacing[3];
  inData->GetSpacing( dSpacing ); // get the spacing

  std::cout << "Spacing = " << dSpacing[0] << " " << dSpacing[1] << " " << dSpacing[2] << std::endl;

  double dOrigin[3];
  inData->GetOrigin( dOrigin ); // get the origin

  std::cout << "Origin = " << dOrigin[0] << " " << dOrigin[1] << " " << dOrigin[2] << std::endl;

  // get memory for the topologcially corrected volume

  if ( pCorrectedImageData != NULL ) {pCorrectedImageData->Delete(); }
  pCorrectedImageData = vtkImageData::New();

  int iExtent[6];
  inData->GetExtent( iExtent );
  pCorrectedImageData->SetExtent( iExtent );
  pCorrectedImageData->SetSpacing( dSpacing );

  // set up and run the genus 0 code

  genus0parameters g0[1]; /* need an instance of genus0parameters */
  genus0init(g0);         /* initialize the instance, set default parameters */

  int            totlen;
  unsigned short *input;
  /* set g0->dims[0..2] and allocate memory */
  totlen = 1;
  for ( int iI = 0; iI < 3; iI++ )
    {
    totlen *= ( g0->dims[iI] = dims[iI] );
    }

  // allocate the temporary memory
  input = (unsigned short *)calloc( totlen, sizeof( unsigned short ) );

  void *ptr = inData->GetScalarPointer();

  switch ( inData->GetScalarType() )
    {
    vtkTemplateMacro(
      vtkGetUnsignedShortData(this, inData, static_cast<VTK_TT *>( ptr ), input, totlen );
      );
    default:
      vtkErrorMacro(<< "Unknown output ScalarType");
      return;
    }

  float ijk2ras[16];
  for ( int iI = 0; iI < 16; iI++ )
    {
    ijk2ras[iI] = 0.0;
    }
  ijk2ras[0] = (float)dSpacing[0];
  ijk2ras[5] = (float)dSpacing[1];
  ijk2ras[10] = (float)dSpacing[2];
  ijk2ras[15] = 1;

  /*float ijk2ras[16];

  for ( int iI=0; iI<4; iI++ ) {
    for ( int iJ=0; iJ<4; iJ++ ) {

      int iC = iI*4+iJ;
      ijk2ras[iC] = matIJKtoRAS->GetElement(iI,iJ);

    }
    }*/

  // and set the parameters to run the code

  // if ( Verbose ) {
  if ( 1 )
    {
    std::cout << "Using verbose mode:" << std::endl << std::endl;
    std::cout << "BiggestComponent = " << BiggestComponent << std::endl;
    std::cout << "ConnectedComponent = " << ConnectedComponent << std::endl;
    std::cout << "CutLoops = " << CutLoops << std::endl;
    std::cout << "ComputeSurface = " << ComputeSurface << std::endl;
    std::cout << "iConnectivity = " << iConnectivity << std::endl;
    }

  g0->connectivity = iConnectivity;
  g0->connected_component = ConnectedComponent;
  g0->input = input;
  g0->value = 1;

  int DesiredCutLoopsValue = CutLoops;

  // There are a number of options which can be computed and are
  // currently supported
  //
  // 6 connectivity: volume and surface computations
  //                 in case of surface computation we compute on the
  //                 inverse volume for the cut-loops option
  //                 Volume computation results in 6-connected
  //                 with possibly multiple connected components
  //                 (requires (if desired) the extraction of the
  //                 largest CC)
  //
  // 18 connectivity: only surfaces are supported currently
  //

  if ( CutLoops && iConnectivity == 6 && !ComputeSurface )
    {
    g0->value = 0;
    DesiredCutLoopsValue = 0;
    SetAltValue( 0 );
    }
  else if ( !CutLoops && iConnectivity == 6 && !ComputeSurface )
    {
    g0->value = 1;
    DesiredCutLoopsValue = 0;
    SetAltValue( 1 );
    }
  else
    {
    g0->value = 1;
    SetAltValue( 1 );
    }

  g0->alt_value = 1;
  g0->contour_value = 1;
  g0->alt_contour_value = 1;

  g0->cut_loops = DesiredCutLoopsValue;

  g0->alt_value = altValue;
  g0->contour_value = altValue;
  g0->alt_contour_value = altValue;

  std::cout << "g0->cut_loops = " << g0->cut_loops << std::endl;
  std::cout << "g0->value = " << g0->value << std::endl;
  std::cout << "g0->alt_value = " << g0->value << std::endl;
  std::cout << "g0->contour_value = " << g0->value << std::endl;
  std::cout << "g0->alt_contour_value = " << g0->value << std::endl;

  g0->any_genus = 0;
  g0->biggest_component = BiggestComponent;
  g0->pad[0] = g0->pad[1] = g0->pad[2] = 2;
  g0->ijk2ras = ijk2ras;
  g0->verbose = Verbose;
  g0->extraijkscale[2] = 1;

  /* call the function! */
  if ( genus0(g0) )
    {
    std::cerr << "Error when executing genus0." << std::endl;
    }

  // now write everything to the polygonal output structures

  int estimatedPoints = 10000;
  int estimatedTriangles = 5000;

  vtkPoints *_Points = vtkPoints::New();
  _Points->Allocate(estimatedPoints, estimatedPoints / 2);
  vtkCellArray *_Triangles = vtkCellArray::New();
  _Triangles->Allocate(estimatedTriangles, estimatedTriangles / 2);

  // write it out

  for ( int iI = 0; iI < g0->vert_count; iI++ )
    {
    double cCoors[3];
    cCoors[0] = g0->vertices[iI];                      // -dOrigin[0];
    cCoors[1] = g0->vertices[iI + g0->vert_count];     // -dOrigin[1];
    cCoors[2] = g0->vertices[iI + 2 * g0->vert_count]; // +dOrigin[2];
    _Points->InsertNextPoint( cCoors );
    }

  for ( int iI = 0; iI < g0->tri_count; iI++ )
    {
    _Triangles->InsertNextCell( 3 );
    _Triangles->InsertCellPoint( g0->triangles[iI] );
    _Triangles->InsertCellPoint( g0->triangles[iI + g0->tri_count] );
    _Triangles->InsertCellPoint( g0->triangles[iI + 2 * g0->tri_count] );
    }

  // and associate it with the output

  outData->SetPoints(_Points);
  _Points->Delete();
  _Points = NULL;
  outData->SetPolys(_Triangles);
  _Triangles->Delete();
  _Triangles = NULL;

  // now export the topologically corrected image
  //g0->output
  void *ptrCI = pCorrectedImageData->GetScalarPointer();

  switch ( pCorrectedImageData->GetScalarType() )
    {
    vtkTemplateMacro(
      vtkSetUnsignedShortData(this, pCorrectedImageData, static_cast<VTK_TT *>( ptrCI ), g0->output, totlen )
      );
    default:
      vtkErrorMacro(<< "Unknown output ScalarType");
      return;
    }

  // determine the number of connected compontens by counting the number of
  // different labels in the output

  std::map<int, int> mapLabels;

  for ( int iI = 0; iI < totlen; iI++ )
    {
    mapLabels[( g0->output )[iI]] = 1;
    }

  this->iConnectedComponents = mapLabels.size() - 1;  // remove one because of
                                                      // the zero label

  // and free everything that is left
  genus0destruct(g0);

  // deallocate the temporary memory
  free(input);
}

// ----------------------------------------------------------------------------
int vtkImageGenus0MarchingCubes::FillInputPortInformation(int, vtkInformation *info)
{
  info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkImageData");
  return 1;
}

// ----------------------------------------------------------------------------
void vtkImageGenus0MarchingCubes::PrintSelf(ostream & os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);
}
