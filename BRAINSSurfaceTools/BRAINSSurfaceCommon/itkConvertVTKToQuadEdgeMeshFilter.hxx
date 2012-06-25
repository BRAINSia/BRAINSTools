/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkConvertVTKToQuadEdgeMeshFilter.txx,v $
  Language:  C++
  Date:      $Date: 2009-06-02 12:48:35 $
  Version:   $Revision: 1.16 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkConvertVTKToQuadEdgeMeshFilter_hxx
#define __itkConvertVTKToQuadEdgeMeshFilter_hxx

#include "itkConvertVTKToQuadEdgeMeshFilter.h"

namespace itk
{
//
// Constructor
//
template <class TOutputMesh>
ConvertVTKToQuadEdgeMeshFilter<TOutputMesh>
::ConvertVTKToQuadEdgeMeshFilter()
{
  //
  // Create the output
  //
  typename TOutputMesh::Pointer output = TOutputMesh::New();
  this->ProcessObject::SetNumberOfRequiredOutputs(1);
  this->ProcessObject::SetNthOutput(0, output.GetPointer() );
}

//
// Set the polydata
//
template <class TOutputMesh>
void
ConvertVTKToQuadEdgeMeshFilter<TOutputMesh>
::SetPolyData(vtkPolyData * polyData)
{
  m_inputPolyData = polyData;
}

//
// Get the polydata
//
template <class TOutputMesh>
vtkPolyData *
ConvertVTKToQuadEdgeMeshFilter<TOutputMesh>
::GetPolyData()
{
  return m_inputPolyData;
}

template <class TOutputMesh>
void
ConvertVTKToQuadEdgeMeshFilter<TOutputMesh>
::GenerateData()
{
  typename OutputMeshType::Pointer outputMesh = this->GetOutput();
  outputMesh->SetCellsAllocationMethod( OutputMeshType::CellsAllocatedDynamicallyCellByCell );

  //
  // Load the point coordinates into the itk::Mesh
  //
  int numberOfPoints = m_inputPolyData->GetNumberOfPoints();
  outputMesh->GetPoints()->Reserve( numberOfPoints );

  PointType point;
  double    vtkPoint[3];
  for( int i = 0; i < numberOfPoints; i++ )
    {
    m_inputPolyData->GetPoint(i, vtkPoint);
    point[0] = vtkPoint[0];
    point[1] = vtkPoint[1];
    point[2] = vtkPoint[2];

    outputMesh->SetPoint( i, point );
    }

  //
  // Load the polygons into the itk::Mesh
  //
  int numberOfCells = m_inputPolyData->GetNumberOfCells();

  vtkIdList *pointList = vtkIdList::New();
  for( int i = 0; i < numberOfCells; i++ )
    {
    CellAutoPointer cell;

    TriangleCellType * triangleCell = new TriangleCellType;
    int                numberOfCellPoints = 3;

    m_inputPolyData->GetCellPoints(i, pointList);
    for( int k = 0; k < numberOfCellPoints; k++ )
      {
      int pointId = pointList->GetId(k);
      triangleCell->SetPointId( k, pointId );
      }

    cell.TakeOwnership( triangleCell );
    outputMesh->SetCell( i, cell );
    }
  pointList->Delete();

  //
  // Load the PointData into the itk::Mesh
  //
  vtkPointData * inputPointData = m_inputPolyData->GetPointData();
  if( inputPointData != NULL )
    {
    vtkDataArray * dataArray = m_inputPolyData->GetPointData()->GetScalars();

    if( dataArray != NULL )
      {
      typedef typename OutputMeshType::PointDataContainer PointDataContainer;

      outputMesh->SetPointData( PointDataContainer::New() );
      outputMesh->GetPointData()->Reserve( numberOfPoints );

      // Read the scalar data
      typename OutputMeshType::PixelType pointData;
      for( int i = 0; i < numberOfPoints; i++ )
        {
        pointData = static_cast<PixelType>( dataArray->GetTuple1(i) );
        outputMesh->SetPointData( i, pointData );
        }
      }
    }
}
} // end of namespace itk

#endif
