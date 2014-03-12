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

 Program:   GTRACT (Guided Tensor Restore Anatomical Connectivity Tractography)
 Module:    $RCSfile: $
 Language:  C++
 Date:      $Date: 2006/03/29 14:53:40 $
 Version:   $Revision: 1.9 $

   Copyright (c) University of Iowa Department of Radiology. All rights reserved.
   See GTRACT-Copyright.txt or http://mri.radiology.uiowa.edu/copyright/GTRACT-Copyright.txt
   for details.

      This software is distributed WITHOUT ANY WARRANTY; without even
      the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
      PURPOSE.  See the above copyright notices for more information.

=========================================================================*/

#include <iostream>
#include <fstream>

#include <vtkPoints.h>
#include <vtkFloatArray.h>
#include <vtkPointData.h>
#include <vtkLookupTable.h>
#include <vtkCellArray.h>
#include <vtkIdList.h>
#include <vtkPolyData.h>
#include <vtkPolyDataWriter.h>
#include <vtkPolyDataReader.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkXMLPolyDataReader.h>
// #include <vtkXMLImageDataWriter.h>
#include <vtkAppendPolyData.h>
#include <vtkSplineFilter.h>
#include "vnl/vnl_math.h"

// ///////////// VTK Version Compatibility   //////////////////////////////
#ifndef vtkFloatingPointType
#define vtkFloatingPointType vtkFloatingPointType
typedef float vtkFloatingPointType;
#endif
// ////////////////////////////////////////////////////////////////////////

#include "compareTractInclusionCLP.h"
#include "BRAINSThreadControl.h"
#include <BRAINSCommonLib.h>

double PairOffFibers(vtkPolyData *resampledTestFibers, vtkPolyData *resampledStandardFibers, int numberOfPoints)
{
  double maxDistances = 0.0;

  for( int j = 0; j < resampledTestFibers->GetNumberOfCells(); j++ )
    {
    if( resampledTestFibers->GetCellType(j) == VTK_POLY_LINE )
      {
      vtkIdList *testPointList = vtkIdList::New();
      resampledTestFibers->GetCellPoints(j, testPointList);
      double minDist = 1E200;
      int    closestK = -1;
      for( int k = 0; k < resampledStandardFibers->GetNumberOfCells(); k++ )
        {
        if( resampledStandardFibers->GetCellType(k) == VTK_POLY_LINE )
          {
          vtkIdList *standardPointList = vtkIdList::New();
          resampledStandardFibers->GetCellPoints(k, standardPointList);
          double sumDist = 0.0;
          for( int i = 0; i < numberOfPoints; i++ )
            {
            vtkFloatingPointType testPoint[3];
            resampledTestFibers->GetPoint(testPointList->GetId(i), testPoint);

            vtkFloatingPointType standardPoint[3];
            resampledStandardFibers->GetPoint(standardPointList->GetId(i), standardPoint);

            double sumSquares = 0.0;
            for( int p = 0; p < 3; p++ )
              {
              double edge = testPoint[p] - standardPoint[p];
              sumSquares += edge * edge;
              }
            sumDist += vcl_sqrt(sumSquares);
            }

          double dist = ( sumDist / numberOfPoints );
          if( dist < minDist )
            {
            minDist = dist;
            closestK = k;
            }
          }
        }

      std::cout << "Pairing test fiber " << j << " with standard fiber " << closestK << " at distance " << minDist
                << std::endl;
      if( maxDistances < minDist )
        {
        maxDistances = minDist;
        }
      }
    }
  return maxDistances;
}

int main(int argc, char * argv[])
{
  PARSE_ARGS;
  BRAINSRegisterAlternateIO();
  const BRAINSUtils::StackPushITKDefaultNumberOfThreads TempDefaultNumberOfThreadsHolder(numberOfThreads);

  const bool debug = true;
  if( debug )
    {
    std::cout << "Test Fiber Tract: " <<  testFiber << std::endl;
    std::cout << "Standard Fiber Tract: " <<  standardFiber << std::endl;
    std::cout << "Use XML PolyData Files: " << writeXMLPolyDataFile  << std::endl;
    std::cout << "Single Fiber Tract Closeness: " <<  closeness << std::endl;
    std::cout << "Test For Tract Bijection: " << testForBijection  << std::endl;
    std::cout << "Test For Tract Cardinality Agreement: " << testForFiberCardinality  << std::endl;
    std::cout << "Number Of Guide Fiber Points: " <<  numberOfPoints << std::endl;
    }

  vtkPolyData *testFiberTract;
  if( writeXMLPolyDataFile )
    {
    vtkXMLPolyDataReader *tractReader = vtkXMLPolyDataReader::New();
    tractReader->SetFileName( testFiber.c_str() );
    tractReader->Update();
    testFiberTract = tractReader->GetOutput();
    }
  else
    {
    vtkPolyDataReader *tractReader = vtkPolyDataReader::New();
    tractReader->SetFileName( testFiber.c_str() );
    tractReader->Update();
    testFiberTract = tractReader->GetOutput();
    }

  vtkPolyData *standardFiberTract;
  if( writeXMLPolyDataFile )
    {
    vtkXMLPolyDataReader *tractReader = vtkXMLPolyDataReader::New();
    tractReader->SetFileName( standardFiber.c_str() );
    tractReader->Update();
    standardFiberTract = tractReader->GetOutput();
    }
  else
    {
    vtkPolyDataReader *tractReader = vtkPolyDataReader::New();
    tractReader->SetFileName( standardFiber.c_str() );
    tractReader->Update();
    standardFiberTract = tractReader->GetOutput();
    }

  int numberOfTestFibers = testFiberTract->GetNumberOfCells();
  int numberOfStandardFibers = standardFiberTract->GetNumberOfCells();
  if( testForFiberCardinality )
    {
    if( numberOfTestFibers != numberOfStandardFibers )
      {
      std::cout << "Number of Test Fibers in " << testFiber << " was " << numberOfTestFibers << std::endl;
      std::cout << "Number of Standard Fibers in " << standardFiber << " was " << numberOfStandardFibers << std::endl;
      std::cout << "TractInclusion test halting with error status 2 indicating cardinality agreement error."
                << std::endl;
      return EXIT_FAILURE;
      }
    }

  vtkSplineFilter *testSpline = vtkSplineFilter::New();
  testSpline->SetInput( testFiberTract );
  testSpline->SetSubdivideToSpecified();
  testSpline->SetNumberOfSubdivisions( numberOfPoints );
  testSpline->Update();

  vtkPolyData *resampledTestFibers = testSpline->GetOutput();

  vtkSplineFilter *standardSpline = vtkSplineFilter::New();
  standardSpline->SetInput( standardFiberTract );
  standardSpline->SetSubdivideToSpecified();
  standardSpline->SetNumberOfSubdivisions( numberOfPoints );
  standardSpline->Update();

  vtkPolyData *resampledStandardFibers = standardSpline->GetOutput();

  double maxDistances = PairOffFibers(resampledTestFibers, resampledStandardFibers, numberOfPoints);
  std::cout << "Maximum distance to standard fibers from all test fibers was " << maxDistances
            << " which is required to be <= " << closeness << std::endl;
  if( !( maxDistances <= closeness ) )  // Also fails on nan
    {
    std::cout
      << "TractInclusion test halting with error status 3 indicating Test fibers not included in Standard neighborhood."
      << std::endl;
    return EXIT_FAILURE;
    }

  if( testForBijection )
    {
    maxDistances = PairOffFibers(resampledStandardFibers, resampledTestFibers, numberOfPoints);
    std::cout << "Maximum distance to test fibers from all standard fibers was " << maxDistances
              << " which is required to be <= " << closeness << std::endl;
    if( !( maxDistances <= closeness ) )  // Also fails on nan
      {
      std::cout
        <<
        "TractInclusion test halting witherror status 4 indicating Standard fibers not included in Test neighborhood."
        << std::endl;
      return EXIT_FAILURE;
      }
    }
  return EXIT_SUCCESS;
}
