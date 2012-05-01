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

#include "vtkPolyDataReader.h"
#include "vtkXMLPolyDataReader.h"
#include "vtkPolyData.h"
#include "vtkDoubleArray.h"
#include "vtkFloatArray.h"
#include "vtkIntArray.h"
#include "vtkShortArray.h"
#include "vtkCellData.h"
#include "vtkPointData.h"
#include "vtkFieldData.h"
#include "vtkTriangle.h"
#include "vtkIdList.h"
#include "itksys/SystemTools.hxx"

#include "BRAINSMeasureSurfaceCLP.h"

#include <stdlib.h>
#include <sys/stat.h>
#include <time.h>
#include <algorithm>

template <class T>
const T & min( const T & a, const T & b )
{
  return ( a < b ) ? a : b;     // or: return comp(a,b)?a:b; for the comp
                                // version
}

int main( int argc, char * *argv )
{
  PARSE_ARGS;

  std::cout << "Surface Measurement Parameters" << std::endl;
  std::cout << "------------------------------------------------------" << std::endl;
  std::cout << "\tInput Surface: " << inputSurface << std::endl;
  std::cout << "\tWrite CSV File: " << writeCsvFile << std::endl;
  std::cout << "\tCSF File: " << csvFile << std::endl;
  std::cout << "\tWrite XML File: " << writeXmlFile << std::endl;
  std::cout << "\tXML File: " << xmlFile << std::endl;
  std::cout << "\tLabel Array Name: " << arrayName << std::endl;
  std::cout << "\tLabel(s): " << std::endl;
  for( unsigned int i = 0; i < labels.size(); i++ )
    {
    std::cout << "\t\t" << i << ": " << labels[i] << std::endl;
    }
  std::cout << "------------------------------------------------------" << std::endl;

  // Read Surface
  vtkPolyData *surfaceData = vtkPolyData::New();

  std::string extension = itksys::SystemTools::LowerCase(
      itksys::SystemTools::GetFilenameExtension(inputSurface) );

  if( extension == ".vtk" )
    {
    vtkPolyDataReader *reader = vtkPolyDataReader::New();
    reader->SetFileName( inputSurface.c_str() );
    reader->Update();
    surfaceData->DeepCopy( reader->GetOutput() );
    reader->Delete();
    }
  else
    {
    vtkXMLPolyDataReader *reader = vtkXMLPolyDataReader::New();
    reader->SetFileName( inputSurface.c_str() );
    reader->Update();
    surfaceData->DeepCopy( reader->GetOutput() );
    reader->Delete();
    }

  vtkShortArray *labelData =
    vtkShortArray::SafeDownCast( surfaceData->GetCellData()->GetAbstractArray( arrayName.c_str() ) );

  double range[2];
  labelData->GetRange(range, 0);
  unsigned int arraySize = static_cast<unsigned int>( range[1] );
  if( arraySize != labels.size() )
    {
    std::cerr << "Error: #label names should match the number of labels" << std::endl;
    std::cerr << " #label names: " << labels.size() << std::endl;
    std::cerr << " #labels: " << arraySize << std::endl;
    return EXIT_FAILURE;
    }

  vtkDoubleArray *surfaceArea = vtkDoubleArray::New();
  vtkDoubleArray *sulcalArea = vtkDoubleArray::New();
  vtkDoubleArray *gyralArea = vtkDoubleArray::New();
  vtkDoubleArray *corticalThickness = vtkDoubleArray::New();
  vtkDoubleArray *sulcalThickness = vtkDoubleArray::New();
  vtkDoubleArray *gyralThickness = vtkDoubleArray::New();
  vtkDoubleArray *corticalCurvature = vtkDoubleArray::New();
  vtkDoubleArray *sulcalCurvature = vtkDoubleArray::New();
  vtkDoubleArray *gyralCurvature = vtkDoubleArray::New();
  vtkIntArray *   totalCount = vtkIntArray::New();
  vtkIntArray *   gyralCount = vtkIntArray::New();
  vtkIntArray *   sulcalCount = vtkIntArray::New();

  // Add one to the size to match labels with the index
  surfaceArea->SetNumberOfValues(arraySize + 1);
  sulcalArea->SetNumberOfValues(arraySize + 1);
  gyralArea->SetNumberOfValues(arraySize + 1);
  corticalThickness->SetNumberOfValues(arraySize + 1);
  sulcalThickness->SetNumberOfValues(arraySize + 1);
  gyralThickness->SetNumberOfValues(arraySize + 1);
  corticalCurvature->SetNumberOfValues(arraySize + 1);
  sulcalCurvature->SetNumberOfValues(arraySize + 1);
  gyralCurvature->SetNumberOfValues(arraySize + 1);
  totalCount->SetNumberOfValues(arraySize + 1);
  gyralCount->SetNumberOfValues(arraySize + 1);
  sulcalCount->SetNumberOfValues(arraySize + 1);
  for( unsigned int i = 0; i < (arraySize + 1); i++ )
    {
    surfaceArea->SetValue(i, 0.0);
    sulcalArea->SetValue(i, 0.0);
    gyralArea->SetValue(i, 0.0);
    corticalThickness->SetValue(i, 0.0);
    sulcalThickness->SetValue(i, 0.0);
    gyralThickness->SetValue(i, 0.0);
    corticalCurvature->SetValue(i, 0.0);
    sulcalCurvature->SetValue(i, 0.0);
    gyralCurvature->SetValue(i, 0.0);
    totalCount->SetValue(i, 0);
    gyralCount->SetValue(i, 0);
    sulcalCount->SetValue(i, 0);
    }

  vtkIdList *    cellPtIds = vtkIdList::New();
  vtkFloatArray *curvature = vtkFloatArray::SafeDownCast( surfaceData->GetPointData()->GetAbstractArray(
                                                            "Mean_Curvature") );
  if( curvature == NULL )
    {
    std::cout << "curvature is empty" << std::endl;
    return EXIT_FAILURE;
    }
  vtkFloatArray *thickness =
    vtkFloatArray::SafeDownCast( surfaceData->GetPointData()->GetAbstractArray("corticalThickness") );
  if( thickness == NULL )
    {
    std::cout << "thickness is empty" << std::endl;
    return EXIT_FAILURE;
    }
  for( int i = 0; i < labelData->GetNumberOfTuples(); i++ )
    {
    int index = labelData->GetValue(i);
    surfaceData->GetCellPoints(i, cellPtIds);
    double point1[3], point2[3], point3[3];
    surfaceData->GetPoint(cellPtIds->GetId(0), point1);
    surfaceData->GetPoint(cellPtIds->GetId(1), point2);
    surfaceData->GetPoint(cellPtIds->GetId(2), point3);
    double cellCurvature[3];
    double cellThickness[3];
    double triangleArea = vtkTriangle::TriangleArea(point1, point2, point3);
    for( int j = 0; j < 3; j++ )
      {
      cellCurvature[j] = curvature->GetValue( cellPtIds->GetId(j) );
      cellThickness[j] = static_cast<double>( thickness->GetValue( cellPtIds->GetId(j) ) );
      }
    double currentThickness = ( cellThickness[0] + cellThickness[1] + cellThickness[2] ) / 3.0;
    double currentCurvature = ( cellCurvature[0] + cellCurvature[1] + cellCurvature[2] ) / 3.0;

    // Curvature and Thickness are point based.
    totalCount->SetValue(index, totalCount->GetValue(index) + 1);
    corticalThickness->SetValue(index, corticalThickness->GetValue(index) + currentThickness * triangleArea);
    surfaceArea->SetValue(index, surfaceArea->GetValue(index) + triangleArea);
    corticalCurvature->SetValue(index, corticalCurvature->GetValue(index) + currentCurvature * triangleArea);

    if( currentCurvature >= 0.0 )
      {
      gyralThickness->SetValue(index, gyralThickness->GetValue(index) + currentThickness * triangleArea);
      gyralArea->SetValue(index, gyralArea->GetValue(index) + triangleArea);
      gyralCurvature->SetValue(index, gyralCurvature->GetValue(index) + currentCurvature * triangleArea);
      gyralCount->SetValue(index, gyralCount->GetValue(index) + 1);
      }
    else
      {
      sulcalThickness->SetValue(index, sulcalThickness->GetValue(index) + currentThickness * triangleArea);
      sulcalArea->SetValue(index, sulcalArea->GetValue(index) + triangleArea);
      sulcalCurvature->SetValue(index, sulcalCurvature->GetValue(index) + currentCurvature * triangleArea);
      sulcalCount->SetValue(index, sulcalCount->GetValue(index) + 1);
      }
    }
  for( unsigned int i = 0; i < (arraySize + 1); i++ )
    {
    if( surfaceArea->GetValue(i) != 0.0 )
      {
      corticalThickness->SetValue( i, corticalThickness->GetValue(i) / surfaceArea->GetValue(i) );
      corticalCurvature->SetValue( i, corticalCurvature->GetValue(i) / surfaceArea->GetValue(i) );
      }
    else
      {
      corticalThickness->SetValue(i, 0.0);
      corticalCurvature->SetValue(i, 0.0);
      }

    if( gyralArea->GetValue(i) != 0.0 )
      {
      gyralThickness->SetValue( i, gyralThickness->GetValue(i) / gyralArea->GetValue(i) );
      gyralCurvature->SetValue( i, gyralCurvature->GetValue(i) / gyralArea->GetValue(i) );
      }
    else
      {
      gyralThickness->SetValue(i, 0.0);
      gyralCurvature->SetValue(i, 0.0);
      }

    if( sulcalArea->GetValue(i) != 0.0 )
      {
      sulcalThickness->SetValue( i, sulcalThickness->GetValue(i) / sulcalArea->GetValue(i) );
      sulcalCurvature->SetValue( i, sulcalCurvature->GetValue(i) / sulcalArea->GetValue(i) );
      }
    else
      {
      sulcalThickness->SetValue(i, 0.0);
      sulcalCurvature->SetValue(i, 0.0);
      }
    }

  std::string measurementTime = itksys::SystemTools::GetCurrentDateTime("%Y-%m-%d_%H:%M.%S");
  std::string user = itksys::SystemTools::GetEnv("USER");

  if( writeCsvFile )
    {
    bool fileExists = itksys::SystemTools::FileExists( csvFile.c_str() );

    std::ofstream csvfile;
    csvfile.open(csvFile.c_str(), std::ios::app);

    if( fileExists == false )
      {
      csvfile << "Subject-ID, Scan-ID, Date, User, Surface, Region";
      csvfile << ", Total Area, Total Mean Thickness, Total Sum Curvature";
      csvfile << ", Gyral Area, Gyral Mean Thickness, Gyral Mean Curvature";
      csvfile << ", Sulcal Area, Sulcal Mean Thickness, Sulcal Mean Curvature" << std::endl;
      }
    for( unsigned int i = 0; i < (arraySize + 1); i++ )
      {
      if( i == 0 )
        {
        csvfile << subjectId << ", " << scanId << ", " << measurementTime << ", " << user << ", " << inputSurface
                << ", " << "Unlabled";
        }
      else
        {
        csvfile << subjectId << ", " << scanId << ", " << measurementTime << ", " << user << ", " << inputSurface
                << ", " << labels[i - 1];
        }

      csvfile << ", " << surfaceArea->GetValue(i);
      csvfile << ", " << corticalThickness->GetValue(i);
      csvfile << ", " << corticalCurvature->GetValue(i);

      csvfile << ", " << gyralArea->GetValue(i);
      csvfile << ", " << gyralThickness->GetValue(i);
      csvfile << ", " << gyralCurvature->GetValue(i);

      csvfile << ", " << sulcalArea->GetValue(i);
      csvfile << ", " << sulcalThickness->GetValue(i);
      csvfile << ", " << sulcalCurvature->GetValue(i) << std::endl;
      }

    csvfile.close();
    }

  if( writeXmlFile )
    {
    std::cerr << "Error: XML file writing is not yet implemented..." << std::endl;
    }

  int result = EXIT_SUCCESS;

  if( testDepth )
    {
    for( unsigned int i = 1; i < (arraySize + 1); i++ )
      {
      if( ( ( fabs(totalDepthResults[i - 1]) - fabs( corticalThickness->GetValue(i) ) )
            / fabs( corticalThickness->GetValue(i) ) ) >= 0.001 )
        {
        std::cerr << "Error: Expected total thickness (" << corticalThickness->GetValue(i) << ") measured (";
        std::cerr <<  totalDepthResults[i - 1] << " )" << std::endl;
        result = EXIT_FAILURE;
        }
      if( ( ( fabs(gyralDepthResults[i
                                     - 1])
              - fabs( gyralThickness->GetValue(i) ) ) / fabs( gyralThickness->GetValue(i) ) ) >= 0.001 )
        {
        std::cerr << "Error: Expected gyral thickness (" << gyralThickness->GetValue(i) << ") measured (";
        std::cerr <<  gyralDepthResults[i - 1] << " )" << std::endl;
        result = EXIT_FAILURE;
        }
      if( ( ( fabs(sulcalDepthResults[i - 1]) - fabs( sulcalThickness->GetValue(i) ) )
            / fabs( sulcalThickness->GetValue(i) ) ) >= 0.001 )
        {
        std::cerr << "Error: Expected sulcal thickness (" << sulcalThickness->GetValue(i) << ") measured (";
        std::cerr <<  sulcalDepthResults[i - 1] << " )" << std::endl;
        result = EXIT_FAILURE;
        }
      }
    }

  if( testArea )
    {
    for( unsigned int i = 1; i < (arraySize + 1); i++ )
      {
      if( ( (fabs(totalAreaResults[i - 1]) - fabs( surfaceArea->GetValue(i) ) ) / fabs( surfaceArea->GetValue(i) ) ) >=
          0.001 )
        {
        std::cerr << "Error: Expected total area (" << surfaceArea->GetValue(i) << ") measured (";
        std::cerr <<  totalDepthResults[i - 1] << " )" << std::endl;
        result = EXIT_FAILURE;
        }
      if( ( ( fabs(gyralAreaResults[i - 1]) - fabs( gyralArea->GetValue(i) ) ) / fabs( gyralArea->GetValue(i) ) ) >=
          0.001 )
        {
        std::cerr << "Error: Expected gyral area (" << gyralArea->GetValue(i) << ") measured (";
        std::cerr <<  gyralDepthResults[i - 1] << " )" << std::endl;
        result = EXIT_FAILURE;
        }
      if( ( ( fabs(sulcalAreaResults[i - 1]) - fabs( sulcalArea->GetValue(i) ) ) / fabs( sulcalArea->GetValue(i) ) ) >=
          0.001 )
        {
        std::cerr << "Error: Expected sulcal area (" << sulcalArea->GetValue(i) << ") measured (";
        std::cerr <<  sulcalDepthResults[i - 1] << " )" << std::endl;
        result = EXIT_FAILURE;
        }
      }
    }

  if( testCurvature )
    {
    for( unsigned int i = 1; i < (arraySize + 1); i++ )
      {
      if( ( ( fabs(totalCurvatureResults[i - 1]) - fabs( corticalCurvature->GetValue(i) ) )
            / fabs( corticalCurvature->GetValue(i) ) ) >= 0.001 )
        {
        std::cerr << "Error: Expected total curvature (" << corticalCurvature->GetValue(i) << ") measured (";
        std::cerr <<  totalCurvatureResults[i - 1] << " )" << std::endl;
        result = EXIT_FAILURE;
        }
      if( ( ( fabs(gyralCurvatureResults[i - 1]) - fabs( gyralCurvature->GetValue(i) ) )
            / fabs( gyralCurvature->GetValue(i) ) ) >= 0.001 )
        {
        std::cerr << "Error: Expected gyral curvature (" << gyralCurvature->GetValue(i) << ") measured (";
        std::cerr <<  gyralCurvatureResults[i - 1] << " )" << std::endl;
        result = EXIT_FAILURE;
        }
      if( ( ( fabs(sulcalCurvatureResults[i - 1]) - fabs( sulcalCurvature->GetValue(i) ) )
            / fabs( sulcalCurvature->GetValue(i) ) ) >= 0.001 )
        {
        std::cerr << "Error: Expected sulcal curvature (" << sulcalCurvature->GetValue(i) << ") measured (";
        std::cerr <<  sulcalCurvatureResults[i - 1] << " )" << std::endl;
        result = EXIT_FAILURE;
        }
      }
    }

  surfaceArea->Delete();
  sulcalArea->Delete();
  gyralArea->Delete();
  corticalThickness->Delete();
  sulcalThickness->Delete();
  gyralThickness->Delete();
  corticalCurvature->Delete();
  sulcalCurvature->Delete();
  gyralCurvature->Delete();
  totalCount->Delete();
  gyralCount->Delete();
  sulcalCount->Delete();
  surfaceData->Delete();
  cellPtIds->Delete();

  return result;
}
