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
#include <vtkAppendPolyData.h>
#include <vtkSplineFilter.h>
#include <BRAINSCommonLib.h>
#include <vtkVersion.h>

// ////////////////////////////////////////////////////////////////////////

#include "gtractCreateGuideFiberCLP.h"
#include "BRAINSThreadControl.h"
int
main(int argc, char * argv[])
{
  PARSE_ARGS;
  BRAINSRegisterAlternateIO();
  const BRAINSUtils::StackPushITKDefaultNumberOfThreads TempDefaultNumberOfThreadsHolder(numberOfThreads);

  const bool debug = true;
  if (debug)
  {
    std::cout << "=====================================================" << std::endl;
    std::cout << "Input Fiber Tract: " << inputFiber << std::endl;
    std::cout << "Output Guide Fiber Tract: " << outputFiber << std::endl;
    std::cout << "Use XML PolyData Files: " << writeXMLPolyDataFile << std::endl;
    std::cout << "Number Of Guide Fiber Points: " << numberOfPoints << std::endl;
    std::cout << "=====================================================" << std::endl;
  }

  vtkPolyData * fiberTract;
  if (writeXMLPolyDataFile)
  {
    vtkXMLPolyDataReader * tractReader = vtkXMLPolyDataReader::New();
    tractReader->SetFileName(inputFiber.c_str());
    tractReader->Update();
    fiberTract = tractReader->GetOutput();
  }
  else
  {
    vtkPolyDataReader * tractReader = vtkPolyDataReader::New();
    tractReader->SetFileName(inputFiber.c_str());
    tractReader->Update();
    fiberTract = tractReader->GetOutput();
  }

  vtkSplineFilter * spline = vtkSplineFilter::New();
#if (VTK_MAJOR_VERSION < 6)
  spline->SetInput(fiberTract);
#else
  spline->SetInputData(fiberTract);
#endif
  spline->SetSubdivideToSpecified();
  spline->SetNumberOfSubdivisions(numberOfPoints);
  spline->Update();

  vtkPolyData * resampledFibers = spline->GetOutput();
  /* Average */
  vtkPoints * guidePoints = vtkPoints::New();
  for (int i = 0; i < numberOfPoints; i++)
  {
    double avgPoint[3];
    avgPoint[0] = 0;
    avgPoint[1] = 0;
    avgPoint[2] = 0;
    int N = 0;
    for (int j = 0; j < resampledFibers->GetNumberOfCells(); j++)
    {
      if (resampledFibers->GetCellType(j) == VTK_POLY_LINE)
      {
        N++;
        vtkIdList * cellPointList = vtkIdList::New();
        resampledFibers->GetCellPoints(j, cellPointList);
        double currentPoint[3];
        resampledFibers->GetPoint(cellPointList->GetId(i), currentPoint);
        avgPoint[0] += currentPoint[0];
        avgPoint[1] += currentPoint[1];
        avgPoint[2] += currentPoint[2];
      }
    }
    avgPoint[0] /= static_cast<double>(N);
    avgPoint[1] /= static_cast<double>(N);
    avgPoint[2] /= static_cast<double>(N);
    guidePoints->InsertNextPoint(avgPoint);
  }
  vtkCellArray * line = vtkCellArray::New();
  line->InsertNextCell(guidePoints->GetNumberOfPoints());
  for (int i = 0; i < guidePoints->GetNumberOfPoints(); i++)
  {
    line->InsertCellPoint(i);
  }
  vtkPolyData * guideFiber = vtkPolyData::New();
  guideFiber->SetPoints(guidePoints);
  guideFiber->SetLines(line);

  if (writeXMLPolyDataFile)
  {
    vtkXMLPolyDataWriter * tractWriter = vtkXMLPolyDataWriter::New();
    tractWriter->SetFileName(outputFiber.c_str());
#if (VTK_MAJOR_VERSION < 6)
    tractWriter->SetInput(guideFiber);
#else
    tractWriter->SetInputData(guideFiber);
#endif
    tractWriter->Update();
  }
  else
  {
    vtkPolyDataWriter * tractWriter = vtkPolyDataWriter::New();
    tractWriter->SetFileName(outputFiber.c_str());
#if (VTK_MAJOR_VERSION < 6)
    tractWriter->SetInput(guideFiber);
#else
    tractWriter->SetInputData(guideFiber);
#endif
    tractWriter->Update();
  }
  return EXIT_SUCCESS;
}
