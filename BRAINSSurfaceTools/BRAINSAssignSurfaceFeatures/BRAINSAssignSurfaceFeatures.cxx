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
#include "vtkFloatArray.h"
#include "vtkPointData.h"
#include "vtkDelaunay3D.h"
#include "vtkUnstructuredGrid.h"
#include "vtkGeometryFilter.h"
#include "vtkCellLocator.h"
#include "vtkPoints.h"
#include "vtkCurvatures.h"

#include <BRAINSCommonLib.h>
#include "BRAINSAssignSurfaceFeaturesCLP.h"

#include <vtkSmartPointer.h>
#include "vtkMath.h"

int main( int argc, char * *argv )
{
  PARSE_ARGS;
  BRAINSRegisterAlternateIO();
  if( inputSurfaceFile == "" )
    {
    std::cerr << "No input file specified" << std::endl;
    return EXIT_FAILURE;
    }
  if( outputSurfaceFile == "" )
    {
    std::cerr << "No output surface file specified" << std::endl;
    return EXIT_FAILURE;
    }
  if( PC.size() != 3 )
    {
    std::cerr << "ERROR: PC point must have three elements" << std::endl;
    return EXIT_FAILURE;
    }
  std::cout << "Assign Surface Features Parameters" << std::endl;
  std::cout << "------------------------------------------------------" << std::endl;
  std::cout << "\tInput Surface FileName: " << inputSurfaceFile << std::endl;
  if( corticalThickness )
    {
    if( outerSurfaceFile == "" )
      {
      std::cerr << "No outer surface file specified" << std::endl;
      return EXIT_FAILURE;
      }
    std::cout << "\tCalculate cortical thickness " << std::endl;
    std::cout << "\tOut-layer Surface: " << outerSurfaceFile << std::endl;
    std::cout << "\tMaximum thickness value: " << maxThickness << std::endl;
    }
  std::cout << "\tOutput Surface: " << outputSurfaceFile << std::endl;

  if( distanceToPC_AP )
    {
    std::cout << "\tCalculate DistanceToPC in the direction of anterior-to-posterior " << std::endl;
    }
  if( distanceToPC_IS )
    {
    std::cout << "\tCalculate DistanceToPC in the direction of inferior-to-superior " << std::endl;
    }
  if( distanceToPC_AP || distanceToPC_IS )
    {
    std::cout << "\tPC: " << PC[0] << " " << PC[1] << " " << PC[2] << std::endl;
    }
  if( distanceToHull )
    {
    std::cout << "\tCalculate DistanceToHull " << std::endl;
    }
  if( curvature )
    {
    std::cout << "\tCalculate Curvature " << std::endl;
    std::cout << "\tCurvature type: " << curvatureType << std::endl;
    }
  std::cout << "------------------------------------------------------" << std::endl;

  // Read InputFile
  vtkSmartPointer<vtkPolyDataReader> surfaceReader = vtkSmartPointer<vtkPolyDataReader>::New();
  surfaceReader->SetFileName(inputSurfaceFile.c_str() );
  surfaceReader->Update();

  vtkSmartPointer<vtkPolyData> surface = surfaceReader->GetOutput();
  int                          npoints = surface->GetNumberOfPoints();

  // DistanceToPC
  if( distanceToPC_AP || distanceToPC_IS )
    {
    double pcPoint[3];
    pcPoint[0] = PC[0];
    pcPoint[1] = PC[1];
    pcPoint[2] = PC[2];

    // go through the points on surface
    // and calculate the distance
    vtkSmartPointer<vtkFloatArray> distArray_IS = vtkSmartPointer<vtkFloatArray>::New();
    vtkSmartPointer<vtkFloatArray> distArray_AP = vtkSmartPointer<vtkFloatArray>::New();

    double pOnSurface[3];
    for( int i = 0; i < npoints; i++ )
      {
      surface->GetPoint(i, pOnSurface);
      int    dir = 0;
      double distPC = 0.0;
      if( distanceToPC_AP )
        {
        dir = 1;

        distPC = pOnSurface[dir] - pcPoint[dir];
        distArray_AP->InsertValue(i, distPC);
        }

      if( distanceToPC_IS )
        {
        dir = 2;

        distPC = pOnSurface[dir] - pcPoint[dir];
        distArray_IS->InsertValue(i, distPC);
        }
      }

    // add the array to surface
    distArray_AP->SetName("distToPC_AP");
    distArray_IS->SetName("distToPC_IS");

    if( distanceToPC_AP )
      {
      if( surface->GetPointData()->GetScalars() == NULL )
        {
        surface->GetPointData()->SetScalars(distArray_AP);
        }
      else
        {
        surface->GetPointData()->AddArray(distArray_AP);
        }
      }

    if( distanceToPC_IS )
      {
      if( surface->GetPointData()->GetScalars() == NULL )
        {
        surface->GetPointData()->SetScalars(distArray_IS);
        }
      else
        {
        surface->GetPointData()->AddArray(distArray_IS);
        }
      }
    }

  // distance to a convex hull of the surface
  if( distanceToHull )
    {
    // generate a convex hull
    vtkSmartPointer<vtkDelaunay3D> del = vtkSmartPointer<vtkDelaunay3D>::New();
    del->SetInput( surface );

    vtkUnstructuredGrid *CurrentMesh = del->GetOutput();

    vtkSmartPointer<vtkGeometryFilter> extractSurface = vtkSmartPointer<vtkGeometryFilter>::New();
    extractSurface->SetInput( CurrentMesh );
    extractSurface->Update();

    vtkSmartPointer<vtkPolyData> hull = extractSurface->GetOutput();

    // calculate the distance from surface to hull
    // create the cell locator
    vtkSmartPointer<vtkCellLocator> cellLocator = vtkSmartPointer<vtkCellLocator>::New();
    cellLocator->SetDataSet(hull);
    cellLocator->BuildLocator();

    vtkIdType                      cellId; int subId;
    double                         distHull, closepoint[3];
    vtkSmartPointer<vtkFloatArray> depthArray = vtkSmartPointer<vtkFloatArray>::New();
    for( int i = 0; i < npoints; i++ )
      {
      // find closest point
      cellLocator->FindClosestPoint(surface->GetPoint(i),
                                    closepoint, cellId, subId, distHull);
      // convert dist2 to dist
      distHull = sqrt(distHull);
      depthArray->InsertValue(i, distHull);
      }

    // add the array to surface
    depthArray->SetName("distToHull");
    if( surface->GetPointData()->GetScalars() == NULL )
      {
      surface->GetPointData()->SetScalars(depthArray);
      }
    else
      {
      surface->GetPointData()->AddArray(depthArray);
      }
    }

  // cortical thickness
  // double distance from 190 to 130 surfaces
  if( corticalThickness )
    {
    // read in the outer surface
    vtkSmartPointer<vtkPolyDataReader> outerSurfaceReader = vtkSmartPointer<vtkPolyDataReader>::New();
    outerSurfaceReader->SetFileName(outerSurfaceFile.c_str() );
    outerSurfaceReader->Update();

    vtkSmartPointer<vtkPolyData> outerSurface = outerSurfaceReader->GetOutput();

    // calculate the distance from input surface to the outer surface
    // create the cell locator
    vtkSmartPointer<vtkCellLocator> cellLocator_thickness = vtkSmartPointer<vtkCellLocator>::New();
    cellLocator_thickness->SetDataSet(outerSurface);
    cellLocator_thickness->BuildLocator();

    vtkIdType                      cellId; int subId;
    double                         thickness, closepoint[3];
    vtkSmartPointer<vtkFloatArray> thickArray = vtkSmartPointer<vtkFloatArray>::New();
    for( int i = 0; i < npoints; i++ )
      {
      // find closest point
      cellLocator_thickness->FindClosestPoint(surface->GetPoint(i),
                                              closepoint, cellId, subId, thickness);
      // convert dist2 to dist
      thickness = sqrt(thickness);
      if( thickness > maxThickness )
        {
        thickness = maxThickness;
        }
      thickArray->InsertValue(i, thickness);
      }

    // add the array to surface
    thickArray->SetName("corticalThickness");
    if( surface->GetPointData()->GetScalars() == NULL )
      {
      surface->GetPointData()->SetScalars(thickArray);
      }
    else
      {
      surface->GetPointData()->AddArray(thickArray);
      }
    }

  // curvature
  if( curvature )
    {
    vtkSmartPointer<vtkCurvatures> curve = vtkSmartPointer<vtkCurvatures>::New();
    curve->SetInput(surface);
    curve->SetCurvatureTypeToMean();
    if( curvatureType == "Gauss" )
      {
      curve->SetCurvatureTypeToGaussian();
      }
    curve->Update();

    // extract curvature values from curve
    // add it as an array on surface
    vtkSmartPointer<vtkFloatArray> curveArray = vtkSmartPointer<vtkFloatArray>::New();
    for( int i = 0; i < npoints; i++ )
      {
      double curveValue = curve->GetOutput()->GetPointData()->GetScalars()->GetTuple1(i);
      curveArray->InsertValue(i, curveValue);
      }

    std::string arrayName = curvatureType + "_Curvature";

    curveArray->SetName(arrayName.c_str() );
    if( surface->GetPointData()->GetScalars() == NULL )
      {
      surface->GetPointData()->SetScalars(curveArray);
      }
    else
      {
      surface->GetPointData()->AddArray(curveArray);
      }
    }

  // write the output surface
  vtkSmartPointer<vtkPolyDataWriter> surfaceWriter = vtkSmartPointer<vtkPolyDataWriter>::New();
  surfaceWriter->SetInput( surface );
  surfaceWriter->SetFileName( outputSurfaceFile.c_str() );
  surfaceWriter->Update();

  return EXIT_SUCCESS;
}
