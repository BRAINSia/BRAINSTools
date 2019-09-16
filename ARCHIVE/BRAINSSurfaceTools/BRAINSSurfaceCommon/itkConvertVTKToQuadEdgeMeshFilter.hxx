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
#ifndef __itkConvertVTKToQuadEdgeMeshFilter_hxx
#define __itkConvertVTKToQuadEdgeMeshFilter_hxx

#include "itkConvertVTKToQuadEdgeMeshFilter.h"

namespace itk
{
//
// Constructor
//
template <typename TOutputMesh>
ConvertVTKToQuadEdgeMeshFilter<TOutputMesh>::ConvertVTKToQuadEdgeMeshFilter()
{
  //
  // Create the output
  //
  typename TOutputMesh::Pointer output = TOutputMesh::New();
  this->ProcessObject::SetNumberOfRequiredOutputs(1);
  this->ProcessObject::SetNthOutput(0, output.GetPointer());
}

//
// Set the polydata
//
template <typename TOutputMesh>
void
ConvertVTKToQuadEdgeMeshFilter<TOutputMesh>::SetPolyData(vtkPolyData * polyData)
{
  m_inputPolyData = polyData;
}

//
// Get the polydata
//
template <typename TOutputMesh>
vtkPolyData *
ConvertVTKToQuadEdgeMeshFilter<TOutputMesh>::GetPolyData()
{
  return m_inputPolyData;
}

template <typename TOutputMesh>
void
ConvertVTKToQuadEdgeMeshFilter<TOutputMesh>::GenerateData()
{
  typename OutputMeshType::Pointer outputMesh = this->GetOutput();
  outputMesh->SetCellsAllocationMethod(OutputMeshType::CellsAllocatedDynamicallyCellByCell);

  //
  // Load the point coordinates into the itk::Mesh
  //
  int numberOfPoints = m_inputPolyData->GetNumberOfPoints();
  outputMesh->GetPoints()->Reserve(numberOfPoints);

  PointType point;
  double    vtkPoint[3];
  for (int i = 0; i < numberOfPoints; i++)
  {
    m_inputPolyData->GetPoint(i, vtkPoint);
    point[0] = vtkPoint[0];
    point[1] = vtkPoint[1];
    point[2] = vtkPoint[2];

    outputMesh->SetPoint(i, point);
  }

  //
  // Load the polygons into the itk::Mesh
  //
  int numberOfCells = m_inputPolyData->GetNumberOfCells();

  vtkIdList * pointList = vtkIdList::New();
  for (int i = 0; i < numberOfCells; i++)
  {
    CellAutoPointer cell;

    TriangleCellType * triangleCell = new TriangleCellType;
    int                numberOfCellPoints = 3;

    m_inputPolyData->GetCellPoints(i, pointList);
    for (int k = 0; k < numberOfCellPoints; k++)
    {
      int pointId = pointList->GetId(k);
      triangleCell->SetPointId(k, pointId);
    }

    cell.TakeOwnership(triangleCell);
    outputMesh->SetCell(i, cell);
  }
  pointList->Delete();

  //
  // Load the PointData into the itk::Mesh
  //
  vtkPointData * inputPointData = m_inputPolyData->GetPointData();
  if (inputPointData != nullptr)
  {
    vtkDataArray * dataArray = m_inputPolyData->GetPointData()->GetScalars();

    if (dataArray != nullptr)
    {
      using PointDataContainer = typename OutputMeshType::PointDataContainer;

      outputMesh->SetPointData(PointDataContainer::New());
      outputMesh->GetPointData()->Reserve(numberOfPoints);

      // Read the scalar data
      typename OutputMeshType::PixelType pointData;
      for (int i = 0; i < numberOfPoints; i++)
      {
        pointData = static_cast<PixelType>(dataArray->GetTuple1(i));
        outputMesh->SetPointData(i, pointData);
      }
    }
  }
}
} // end of namespace itk

#endif
