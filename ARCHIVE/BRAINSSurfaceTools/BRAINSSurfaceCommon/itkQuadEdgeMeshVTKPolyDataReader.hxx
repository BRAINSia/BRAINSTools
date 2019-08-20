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
#ifndef __itkQuadEdgeMeshVTKPolyDataReader_hxx
#define __itkQuadEdgeMeshVTKPolyDataReader_hxx

#include "itkQuadEdgeMeshVTKPolyDataReader.h"
#include <fstream>
#include <cstdio>
#include <cstring>

namespace itk
{
//
// Constructor
//
template < typename TOutputMesh >
QuadEdgeMeshVTKPolyDataReader< TOutputMesh >::QuadEdgeMeshVTKPolyDataReader()
{
  //
  // Create the output
  //
  typename TOutputMesh::Pointer output = TOutputMesh::New();
  this->ProcessObject::SetNumberOfRequiredOutputs( 1 );
  this->ProcessObject::SetNthOutput( 0, output.GetPointer() );
}

template < typename TOutputMesh >
void
QuadEdgeMeshVTKPolyDataReader< TOutputMesh >::GenerateData()
{
  typename OutputMeshType::Pointer outputMesh = this->GetOutput();

  outputMesh->SetCellsAllocationMethod( OutputMeshType::CellsAllocatedDynamicallyCellByCell );

  if ( m_FileName == "" )
  {
    itkExceptionMacro( "No input FileName" );
  }

  //
  // Read input file
  //
  std::ifstream inputFile( m_FileName.c_str() );

  if ( !inputFile.is_open() )
  {
    itkExceptionMacro( "Unable to open file\n"
                       "inputFilename= "
                       << m_FileName );
  }

  std::string line;

  while ( !inputFile.eof() )
  {
    std::getline( inputFile, line );

    if ( line.find( "POINTS" ) != std::string::npos )
    {
      break;
    }
  }

  itkDebugMacro( "POINTS line" << line );

  std::string pointLine( line, strlen( "POINTS " ), line.length() );
  itkDebugMacro( "pointLine " << pointLine );

  int numberOfPoints = -1;

  if ( sscanf( pointLine.c_str(), "%d", &numberOfPoints ) != 1 )
  {
    itkExceptionMacro( "ERROR: Failed to read numberOfPoints\n"
                       "       pointLine= "
                       << pointLine );
  }

  itkDebugMacro( "numberOfPoints= " << numberOfPoints );

  if ( numberOfPoints < 1 )
  {
    itkExceptionMacro( "numberOfPoints < 1"
                       << "       numberOfPoints= " << numberOfPoints );
  }

  outputMesh->GetPoints()->Reserve( numberOfPoints );

  //
  // Load the point coordinates into the itk::Mesh
  //

  PointType point;
  for ( int i = 0; i < numberOfPoints; i++ )
  {
    inputFile >> point;
    outputMesh->SetPoint( i, point );
  }

  // Continue searching for the POLYGONS line
  while ( !inputFile.eof() && line.find( "POLYGONS" ) == std::string::npos )
  {
    std::getline( inputFile, line );
  }

  itkDebugMacro( "POLYGONS line" << line );

  std::string polygonLine( line, strlen( "POLYGONS " ), line.length() );
  itkDebugMacro( "polygonLine " << polygonLine );

  //
  // Read the number of polygons
  //

  CellIdentifier numberOfPolygons = 0;
  CellIdentifier numberOfIndices = 0;

  if ( sscanf( polygonLine.c_str(), "%lu %lu", &numberOfPolygons, &numberOfIndices ) != 2 )
  {
    itkExceptionMacro( "ERROR: Failed to read numberOfPolygons from subline2"
                       "\npolygonLine= "
                       << polygonLine );
  }

  itkDebugMacro( "numberOfPolygons " << numberOfPolygons );
  itkDebugMacro( "numberOfIndices " << numberOfIndices );

  if ( numberOfPolygons < 1 )
  {
    itkExceptionMacro( "ERROR: numberOfPolygons < 1\nnumberOfPolygons= " << numberOfPolygons );
  }

  if ( numberOfIndices < numberOfPolygons )
  {
    itkExceptionMacro( "ERROR: numberOfIndices < numberOfPolygons\n"
                       << "numberOfIndices= " << numberOfIndices << "\n"
                       << "numberOfPolygons= " << numberOfPolygons );
  }

  //
  // Load the polygons into the itk::Mesh
  //

  PointIdentifier numberOfCellPoints;
  long            ids[3];
  for ( CellIdentifier i = 0; i < numberOfPolygons; i++ )
  {
    if ( inputFile.eof() )
    {
      itkExceptionMacro( "Failed to read " << numberOfPolygons << " polygons before the end of file" );
    }

    std::getline( inputFile, line );

    if ( line.find( "DATA" ) != std::string::npos )
    {
      itkExceptionMacro( "Read keyword DATA" );
    }

    int got;
    if ( ( got = sscanf( line.c_str(), "%lu %ld %ld %ld", &numberOfCellPoints, &ids[0], &ids[1], &ids[2] ) ) != 4 )
    {
      itkExceptionMacro( "Error parsing POLYGON cell. Expected 4 items but got " << got << std::endl
                                                                                 << "Line is: " << line );
    }
    if ( numberOfCellPoints != 3 )
    {
      itkExceptionMacro( "ERROR: numberOfCellPoints != 3\n"
                         << "numberOfCellPoints= " << numberOfCellPoints
                         << "itkQuadEdgeMeshVTKPolyDataReader can only read triangles" );
    }

    if ( static_cast< long >( ids[0] ) < 0 || static_cast< long >( ids[1] ) < 0 || static_cast< long >( ids[2] ) < 0 )
    {
      itkExceptionMacro( "ERROR: Incorrect point ids\n"
                         "ids="
                         << ids[0] << " " << ids[1] << " " << ids[2] );
    }

    if ( static_cast< long >( ids[0] ) >= numberOfPoints || static_cast< long >( ids[1] ) >= numberOfPoints ||
         static_cast< long >( ids[2] ) >= numberOfPoints )
    {
      itkExceptionMacro( "ERROR: Incorrect point ids\n"
                         << "ids=" << ids[0] << " " << ids[1] << " " << ids[2] );
    }

    CellAutoPointer    cell;
    TriangleCellType * triangleCell = new TriangleCellType;
    for ( PointIdentifier k = 0; k < numberOfCellPoints; k++ )
    {
      triangleCell->SetPointId( k, ids[k] );
    }

    cell.TakeOwnership( triangleCell );
    outputMesh->SetCell( i, cell );
  }

  bool foundPointData = false;

  while ( !inputFile.eof() )
  {
    std::getline( inputFile, line );

    if ( line.find( "POINT_DATA" ) != std::string::npos )
    {
      foundPointData = true;
      break;
    }
  }

  if ( foundPointData )
  {
    using PointDataContainer = typename OutputMeshType::PointDataContainer;

    outputMesh->SetPointData( PointDataContainer::New() );
    outputMesh->GetPointData()->Reserve( numberOfPoints );

    itkDebugMacro( "POINT_DATA line" << line );

    // The following line should be SCALARS or VECTORS
    if ( !inputFile.eof() )
    {
      std::getline( inputFile, line );
    }
    else
    {
      itkExceptionMacro( "Unexpected end-of-file while trying to read POINT_DATA." );
    }

    if ( line.find( "SCALARS" ) != std::string::npos )
    {
      // Skip the following line, since it should contain the LOOKUP_TABLE string
      if ( !inputFile.eof() )
      {
        std::getline( inputFile, line );
        if ( line.find( "LOOKUP_TABLE" ) == std::string::npos )
        {
          itkExceptionMacro( "Expected LOOKUP_TABLE after SCALARS but didn't find it" );
        }
      }
      else
      {
        itkExceptionMacro( "Unexpected end-of-file while trying to read POINT_DATA." );
      }

      // Read the scalar data
      typename OutputMeshType::PixelType pointData;
      for ( int pid = 0; pid < numberOfPoints; pid++ )
      {
        if ( inputFile.eof() )
        {
          itkExceptionMacro( "Unexpected end-of-file while trying to read POINT_DATA."
                             << "Failed while trying to reading point data for id: " << pid );
        }
        inputFile >> pointData;
        outputMesh->SetPointData( pid, pointData );
      }
    }

    if ( line.find( "VECTORS" ) != std::string::npos )
    {
      // Read the vector data
      typename OutputMeshType::PixelType pointData;
      for ( int pid = 0; pid < numberOfPoints; pid++ )
      {
        if ( inputFile.eof() )
        {
          itkExceptionMacro( "Unexpected end-of-file while trying to read POINT_DATA."
                             << "Failed while trying to reading point data for id: " << pid );
        }
        inputFile >> pointData;
        outputMesh->SetPointData( pid, pointData );
      }
    }
  }

  inputFile.close();
}

template < typename TOutputMesh >
void
QuadEdgeMeshVTKPolyDataReader< TOutputMesh >::PrintSelf( std::ostream & os, Indent indent ) const
{
  Superclass::PrintSelf( os, indent );

  os << indent << "FileName: " << m_FileName << std::endl;
}
} // end of namespace itk

#endif
