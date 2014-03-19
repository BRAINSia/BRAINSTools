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
 Date:      $Date: 2006/03/29 14:53:40 $
 Version:   $Revision: 1.9 $

   Copyright (c) University of Iowa Department of Radiology. All rights reserved.
   See GTRACT-Copyright.txt or http://mri.radiology.uiowa.edu/copyright/GTRACT-Copyright.txt
   for details.

      This software is distributed WITHOUT ANY WARRANTY; without even
      the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
      PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
// INPUT:
//   1) Label map that contains at least part of the input surface
//   2) Surface
//
// OUTPUT:
//   1) A copy of the original surface, but with each of its cell values changed
// to the
//     most-common value of its adjacent pixels

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkPoint.h"

#include "vtkPolyDataReader.h"
#include "vtkPolyDataWriter.h"
#include "vtkXMLPolyDataReader.h"
#include "vtkXMLPolyDataWriter.h"
#include "vtkPolyData.h"
#include "vtkErrorCode.h"
#include "vtkCell.h"
#include "vtkIdList.h"
#include "vtkPoints.h"
#include "vtkPointData.h"
#include "vtkCellData.h"
#include "vtkFieldData.h"
#include "vtkShortArray.h"
#include "vtkTransform.h"
#include "vtkTransformPolyDataFilter.h"
#include "vtkDataSetWriter.h"

#include <vtksys/SystemTools.hxx>

#include "BRAINSApplySurfaceLabelsCLP.h"
#include <BRAINSCommonLib.h>

int main( int argc, char * *argv )
{
  PARSE_ARGS;
  BRAINSRegisterAlternateIO();

  std::cout << "Surface Label Map Generation Parameters: \n";
  std::cout << "-----------------------------------------------\n";
  std::cout << "\tLabel Map: " << inputLabelMap << std::endl;
  std::cout << "\tInput Surface: " << inputSurface << std::endl;
  std::cout << "\tCell Data Name: " << cellDataName << std::endl;
  std::cout << "\tOutput Surface: " << outputSurface << std::endl;
  std::cout << "-----------------------------------------------\n";

  // Read Label Map
  typedef itk::Image<signed short, 3>        LabelMapType;
  typedef itk::ImageFileReader<LabelMapType> LabelMapReaderType;

  LabelMapReaderType::Pointer labelMapReader = LabelMapReaderType::New();
  labelMapReader->SetFileName( inputLabelMap.c_str() );
  try
    {
    labelMapReader->Update();
    }
  catch( itk::ExceptionObject & exp )
    {
    std::cerr << "Exception caught!\n";
    std::cerr << exp << std::endl;
    }

  // Read Surface
  vtkPolyData *surface;
  std::string  fileExtension = vtksys::SystemTools::GetFilenameLastExtension(inputSurface);
  std::string  extension = vtksys::SystemTools::LowerCase(fileExtension);
  if( extension == ".vtk" )
    {
    vtkPolyDataReader *surfaceReader = vtkPolyDataReader::New();
    surfaceReader->SetFileName( inputSurface.c_str() );
    try
      {
      surfaceReader->Update();
      }
    catch( vtkErrorCode & exp )
      {
      std::cerr << "Exception caught!\n";
      std::cerr << exp.GetStringFromErrorCode( vtkErrorCode::FileNotFoundError) << std::endl;
      }
    surface = surfaceReader->GetOutput();
    }
  else
    {
    vtkXMLPolyDataReader *surfaceReader = vtkXMLPolyDataReader::New();
    surfaceReader->SetFileName( inputSurface.c_str() );
    try
      {
      surfaceReader->Update();
      }
    catch( vtkErrorCode & exp )
      {
      std::cerr << "Exception caught!\n";
      std::cerr << exp.GetStringFromErrorCode( vtkErrorCode::FileNotFoundError) << std::endl;
      }
    surface = surfaceReader->GetOutput();
    }

  // Put surface into RAS orientation
  vtkTransform *rasOrientation = vtkTransform::New();
  const double  orientationMatrix[16] = { -1,  0, 0, 0,
                                          0, -1, 0, 0,
                                          0,  0, 1, 0,
                                          0,  0, 0, 0  };
  rasOrientation->SetMatrix( orientationMatrix );

  vtkTransformPolyDataFilter *resampleSurface = vtkTransformPolyDataFilter::New();
#if (VTK_MAJOR_VERSION < 6)
  resampleSurface->SetInput( surface );
#else
  resampleSurface->SetInputData( surface );
#endif
  resampleSurface->SetTransform( rasOrientation );
  resampleSurface->Update();

  surface = resampleSurface->GetOutput();

  // Loop variables
  int                     numCells = surface->GetNumberOfCells();
  vtkIdList *             pointList = vtkIdList::New();
  double                  pnt1[3];
  double                  pnt2[3];
  double                  pnt3[3];
  LabelMapType::PointType convertedPoint;
  LabelMapType::IndexType pixelIndex;
  short                   labelValues[3] = { 0, 0, 0 };
  short                   finalLabelValue;
  vtkShortArray *         cellDataArray = vtkShortArray::New();
  cellDataArray->SetName( cellDataName.c_str() );
  cellDataArray->SetNumberOfValues( numCells );
  // Iterate over the surface
  for( int i = 0; i < numCells; i++ )
    {
    // At each cell obtain all 3 of its points
    surface->GetCellPoints( i, pointList);

    // Copy the cell's points into 3 double*'s
    surface->GetPoint( pointList->GetId( 0 ), pnt1 );
    surface->GetPoint( pointList->GetId( 1 ), pnt2 );
    surface->GetPoint( pointList->GetId( 2 ), pnt3 );

    // Obtain the label(from the label map) of each of the 3 points

    // Convert the VTK point to an ITK point

    convertedPoint.SetElement( 0, pnt1[0] );
    convertedPoint.SetElement( 1, pnt1[1] );
    convertedPoint.SetElement( 2, pnt1[2] );
    labelMapReader->GetOutput()->TransformPhysicalPointToIndex( convertedPoint, pixelIndex );
    labelValues[0] = labelMapReader->GetOutput()->GetPixel( pixelIndex );
    // std::cout << "Point 1: " << convertedPoint << std::endl;
    // std::cout << "Index 1: " << pixelIndex << std::endl;
    // std::cout << "Value 1: " << labelValues[0] << std::endl;
    convertedPoint.SetElement( 0, pnt2[0] );
    convertedPoint.SetElement( 1, pnt2[1] );
    convertedPoint.SetElement( 2, pnt2[2] );
    labelMapReader->GetOutput()->TransformPhysicalPointToIndex( convertedPoint, pixelIndex );
    labelValues[1] = labelMapReader->GetOutput()->GetPixel( pixelIndex );

    convertedPoint.SetElement( 0, pnt3[0] );
    convertedPoint.SetElement( 1, pnt3[1] );
    convertedPoint.SetElement( 2, pnt3[2] );
    labelMapReader->GetOutput()->TransformPhysicalPointToIndex( convertedPoint, pixelIndex );
    labelValues[2] = labelMapReader->GetOutput()->GetPixel( pixelIndex );
    // std::cout << "Point 3: " << convertedPoint << std::endl;
    // std::cout << "Index 3: " << pixelIndex << std::endl;
    // std::cout << "Value 3: " << labelValues[2] << std::endl;
    // std::cout << "===============================================" <<
    // std::endl;

    // If 2+ points have the same label, assign that value to the label;
    // else assign the cell value to that of the initial point
    if( labelValues[1] == labelValues[2] )
      {
      finalLabelValue = labelValues[1];
      }
    else
      {
      finalLabelValue = labelValues[0];
      }

    // Change the value of the cell to finalLabelValue
    cellDataArray->SetValue( i, finalLabelValue );
    }

  vtkCellData *cellData = surface->GetCellData();
  cellData->AddArray( cellDataArray );
  //  vtkAbstractArray *tmpArray = surface->GetCellData()->GetAbstractArray( cellDataName.c_str() );

  // Put the surface back into its original orientation

  vtkTransformPolyDataFilter *revertedSurface = vtkTransformPolyDataFilter::New();
#if (VTK_MAJOR_VERSION < 6)
  revertedSurface->SetInput( resampleSurface->GetOutput() );
#else
  revertedSurface->SetInputData( resampleSurface->GetOutput() );
#endif
  revertedSurface->SetTransform( rasOrientation->GetInverse() );
  revertedSurface->Update();

  surface = revertedSurface->GetOutput();

  // Write out the new surface
  fileExtension = vtksys::SystemTools::GetFilenameLastExtension(outputSurface);
  extension = vtksys::SystemTools::LowerCase(fileExtension);
  if( extension == ".vtk" )
    {
    vtkPolyDataWriter *surfaceWriter = vtkPolyDataWriter::New();
#if (VTK_MAJOR_VERSION < 6)
    surfaceWriter->SetInput( surface );
#else
    surfaceWriter->SetInputData( surface );
#endif
    surfaceWriter->SetFileName( outputSurface.c_str() );
    surfaceWriter->Update();
    surfaceWriter->Delete();
    }
  else
    {
    vtkXMLPolyDataWriter *surfaceWriter = vtkXMLPolyDataWriter::New();
#if (VTK_MAJOR_VERSION < 6)
    surfaceWriter->SetInput( surface );
#else
    surfaceWriter->SetInputData( surface );
#endif
    surfaceWriter->SetFileName( outputSurface.c_str() );
    surfaceWriter->SetDataModeToAscii();
    surfaceWriter->Update();
    surfaceWriter->Delete();
    }
  // Clean-up Allocated Objects
  cellDataArray->Delete();
  pointList->Delete();

  return EXIT_SUCCESS;
}
