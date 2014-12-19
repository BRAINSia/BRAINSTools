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
#include "vtkSmartPointer.h"
#include "vtkPointData.h"
#include "vtkPolyDataConnectivityIDFilter.h"
#include "vtkIdList.h"
#include "vtkDataArray.h"
#include "vtkMaskLabel.h"
#include "vtkVersion.h"
#include "SurfaceLabelCleanUpCLP.h"

#include "itkMacro.h" //Needed for ITK_NULLPTR

int SurfaceConnectivityCells(vtkSmartPointer<vtkPolyData> mesh);

int SurfaceConnectivityPoints(vtkSmartPointer<vtkPolyData> mesh);

int RemoveIsolatedPoints(vtkSmartPointer<vtkPolyData> mesh);

void FlipSharpTriangles(vtkSmartPointer<vtkPolyData> mesh);

int main( int argc, char * argv[] )
{
  PARSE_ARGS
  std::cout << "---------------------------------------------------" << std::endl;
  std::cout << "Input Surface: " << std::endl;
  std::cout << inputSurfaceFile << std::endl;
  std::cout << "Output Surface: " << std::endl;
  std::cout << outputSurfaceFile << std::endl;
  std::cout << "---------------------------------------------------" << std::endl;

  // Create all of the classes we will need
  vtkSmartPointer<vtkPolyDataReader> reader = vtkSmartPointer<vtkPolyDataReader>::New();
  reader->SetFileName(inputSurfaceFile.c_str() );
  reader->Update();

  vtkSmartPointer<vtkPolyData> surface_in = reader->GetOutput();
  vtkDataArray *               labelArray = surface_in->GetPointData()->GetScalars();
  std::string                  arrayName = labelArray->GetName();
  if( labelArray == ITK_NULLPTR || arrayName != "LabelValue" )
    {
    std::cerr << "There is no label array as scalars on input surface. ";
    std::cerr << "Quit." << std::endl;
    return 1;
    }

  bool change = true;
  int  iter = 0;
  // iterate until no change in either SurfaceConnectivity or RemoveIsolatedPoints
  while( change && iter < 10 )
    {
    change = false;

    // clean up extra connected regions
    int numConnectedCells = SurfaceConnectivityCells(surface_in);
    int numConnectedPoints = SurfaceConnectivityPoints(surface_in);
    int numIsolated = RemoveIsolatedPoints(surface_in);
    FlipSharpTriangles(surface_in);

    if( numConnectedCells > 0 || numConnectedPoints > 0 || numIsolated > 0 )
      {
      change = true;
      std::cout << "iteration: " << iter << std::endl;
      iter += 1;
      }
    }

  vtkSmartPointer<vtkPolyDataWriter> writer = vtkSmartPointer<vtkPolyDataWriter>::New();
#if (VTK_MAJOR_VERSION < 6)
  writer->SetInput(surface_in);
#else
  writer->SetInputData(surface_in);
#endif
  writer->SetFileName(outputSurfaceFile.c_str() );
  writer->Update();

  return 0;
}

// -------------------------------------------------------------------------------------
// The function cleans up connected regions for each label
// Keep the biggest region and remove any other regions
// Return the number of changes happened to the mesh
// only consider connected cells with all points labeled as mask label
int SurfaceConnectivityCells(vtkSmartPointer<vtkPolyData> mesh)
{
  // initialize number of changes made in this function
  int numChanges = 0;

  // ready for GetPointCells
  mesh->BuildLinks();

  vtkDataArray *labelArray = mesh->GetPointData()->GetScalars();

  // label range
  double labelRange[2];
  labelArray->GetRange(labelRange);

  // set up an array to keep frequency of each label
  int                                              nLabels = int(labelRange[1] - labelRange[0] + 1);
  int *                                            neighborLabels = new int[nLabels];
  vtkSmartPointer<vtkMaskLabel>                    mask = vtkSmartPointer<vtkMaskLabel>::New();
  vtkSmartPointer<vtkPolyDataConnectivityIDFilter> connect = vtkSmartPointer<vtkPolyDataConnectivityIDFilter>::New();
  vtkSmartPointer<vtkPolyData>                     island = vtkSmartPointer<vtkPolyData>::New();
  for( int i = 0; i < nLabels; i++ )
    {
    int label_i = labelRange[0] + i;

    // mask out the label regions first
    // mask filter keeps the number of points
#if (VTK_MAJOR_VERSION < 6)
    mask->SetInput(mesh);
#else
    mask->SetInputData(mesh);
#endif
    mask->LabelOnlyOn();
    mask->SetLabel(label_i);
    mask->Update();

    // extract connected regions for the masked regions
#if (VTK_MAJOR_VERSION < 6)
    connect->SetInput(mask->GetOutput() );
#else
    connect->SetInputData(mask->GetOutput() );
#endif
    connect->ColorRegionsOn();
    connect->SetExtractionModeToAllRegions();
    connect->Update();

    int nRegions = connect->GetNumberOfExtractedRegions();
    // std::cout<<"Label: "<<label_i<<" has "<<nRegions<<" regions."<<std::endl;
    int largestId = connect->GetLargestRegionId();
    // look at each small regions (islands)
    for( int j = 0; j < nRegions; j++ )
      {
      if( j != largestId )
        {
        connect->InitializeSpecifiedRegionList();
        connect->AddSpecifiedRegion(j);
        connect->SetExtractionModeToSpecifiedRegions();

        island = connect->GetOutput();
#if (VTK_MAJOR_VERSION < 6)
        island->Update();
#endif
        island->BuildLinks();

        int ncells = island->GetNumberOfCells();
        // look at each point's LabelValue
        // if it has maximum #of not being "the Label",
        // set the non-label to all of the points on the island.
        // clean up neighborLabels first
        for( int ii = 0; ii < nLabels; ii++ )
          {
          neighborLabels[ii] = 0;
          }

        vtkIdType  npts;
        vtkIdType *pts;
        double     pt[3];
        // put borderIds into a list
        // Ids are from original surface
        vtkSmartPointer<vtkIdList> borderIds = vtkSmartPointer<vtkIdList>::New();
        for( int ii = 0; ii < ncells; ii++ )
          {
          island->GetCellPoints(ii, npts, pts);

          // decide if the point is on the border
          // check neighbor cells of each edge
          vtkSmartPointer<vtkIdList> edgeCells0 = vtkSmartPointer<vtkIdList>::New();
          vtkSmartPointer<vtkIdList> edgeCells1 = vtkSmartPointer<vtkIdList>::New();
          vtkSmartPointer<vtkIdList> edgeCells2 = vtkSmartPointer<vtkIdList>::New();
          island->GetCellEdgeNeighbors(ii, pts[0], pts[1], edgeCells0);
          island->GetCellEdgeNeighbors(ii, pts[0], pts[2], edgeCells1);
          island->GetCellEdgeNeighbors(ii, pts[1], pts[2], edgeCells2);
          // if any edge has zero neighbors
          // the edge is on the border

          int pId_orig;
          if( edgeCells0->GetNumberOfIds() == 0 )
            {
            island->GetPoint(pts[0], pt);
            pId_orig = mesh->FindPoint(pt[0], pt[1], pt[2]);
            borderIds->InsertUniqueId(pId_orig);
            island->GetPoint(pts[1], pt);
            pId_orig = mesh->FindPoint(pt[0], pt[1], pt[2]);
            borderIds->InsertUniqueId(pId_orig);
            }
          if( edgeCells1->GetNumberOfIds() == 0 )
            {
            island->GetPoint(pts[0], pt);
            pId_orig = mesh->FindPoint(pt[0], pt[1], pt[2]);
            borderIds->InsertUniqueId(pId_orig);
            island->GetPoint(pts[2], pt);
            pId_orig = mesh->FindPoint(pt[0], pt[1], pt[2]);
            borderIds->InsertUniqueId(pId_orig);
            }
          if( edgeCells2->GetNumberOfIds() == 0 )
            {
            island->GetPoint(pts[1], pt);
            pId_orig = mesh->FindPoint(pt[0], pt[1], pt[2]);
            borderIds->InsertUniqueId(pId_orig);
            island->GetPoint(pts[2], pt);
            pId_orig = mesh->FindPoint(pt[0], pt[1], pt[2]);
            borderIds->InsertUniqueId(pId_orig);
            }
          }

        // find neighbor Ids for the island
        vtkSmartPointer<vtkIdList> neighborIds = vtkSmartPointer<vtkIdList>::New();
        for( int jj = 0; jj < borderIds->GetNumberOfIds(); jj++ )
          {
          // look at the cells of border point on original surface
          int                        pId_orig = borderIds->GetId(jj);
          vtkSmartPointer<vtkIdList> cellIds = vtkSmartPointer<vtkIdList>::New();

          mesh->GetPointCells(pId_orig, cellIds);
          for( int cell_i = 0; cell_i < cellIds->GetNumberOfIds(); cell_i++ )
            {
            mesh->GetCellPoints(cellIds->GetId(cell_i), npts, pts);
            for( int cellp = 0; cellp < npts; cellp++ )
              {
              if( pts[cellp] != pId_orig )
                {
                neighborIds->InsertUniqueId(pts[cellp]);
                }
              }
            }
          }
        // add neighbor labels
        for( int neighbor_i = 0; neighbor_i < neighborIds->GetNumberOfIds(); neighbor_i++ )
          {
          int label_n = labelArray->GetTuple1(neighborIds->GetId(neighbor_i) );
          if( label_n != label_i )
            {
            neighborLabels[label_n - int(labelRange[0])] += 1;
            }
          }

        // find the newLabel to be the one in neighborLabels with
        // maximum frequency
        int newLabel = label_i;
        for( int ii = 0; ii < nLabels; ii++ )
          {
          if( neighborLabels[ii] > neighborLabels[newLabel - int(labelRange[0])] )
            {
            newLabel = ii + labelRange[0];
            }
          }

        if( newLabel != label_i )
          {
          numChanges += 1;
          }
        // if found newLabel, assign it to all of the points of the island
        // on mesh
        for( int ii = 0; ii < ncells; ii++ )
          {
          island->GetCellPoints(ii, npts, pts);
          for( int jj = 0; jj < npts; jj++ )
            {
            // get each point
            island->GetPoint(pts[jj], pt);
            int pid_orig = mesh->FindPoint(pt[0], pt[1], pt[2]);
            labelArray->SetTuple1(pid_orig, newLabel);
            }
          }
        }
      }
    }
  delete [] neighborLabels;
  return numChanges;
}

// -------------------------------------------------------------------------------------
// The function cleans up connected regions for each label
// Keep the biggest region and remove any other regions
// Return the number of changes happened to the mesh
// Consider connected cells with any point labeled as mask label
int SurfaceConnectivityPoints(vtkSmartPointer<vtkPolyData> mesh)
{
  // initialize number of changes made in this function
  int numChanges = 0;

  // ready for GetPointCells
  mesh->BuildLinks();

  vtkDataArray *labelArray = mesh->GetPointData()->GetScalars();

  // label range
  double labelRange[2];
  labelArray->GetRange(labelRange);

  // set up an array to keep frequency of each label
  int                                              nLabels = int(labelRange[1] - labelRange[0] + 1);
  int *                                            neighborLabels = new int[nLabels];
  vtkSmartPointer<vtkMaskLabel>                    mask = vtkSmartPointer<vtkMaskLabel>::New();
  vtkSmartPointer<vtkPolyDataConnectivityIDFilter> connect = vtkSmartPointer<vtkPolyDataConnectivityIDFilter>::New();
  vtkSmartPointer<vtkPolyData>                     island = vtkSmartPointer<vtkPolyData>::New();
  for( int i = 0; i < nLabels; i++ )
    {
    int label_i = labelRange[0] + i;

    // mask out the label regions first
    // mask filter keeps the number of points
#if (VTK_MAJOR_VERSION < 6)
    mask->SetInput(mesh);
#else
    mask->SetInputData(mesh);
#endif
    mask->SetLabel(label_i);
    mask->Update();

    // extract connected regions for the masked regions
#if (VTK_MAJOR_VERSION < 6)
    connect->SetInput(mask->GetOutput() );
#else
    connect->SetInputData(mask->GetOutput() );
#endif
    connect->ColorRegionsOn();
    connect->SetExtractionModeToAllRegions();
    connect->Update();

    int nRegions = connect->GetNumberOfExtractedRegions();
    // std::cout<<"Label: "<<label_i<<" has "<<nRegions<<" regions."<<std::endl;
    int largestId = connect->GetLargestRegionId();
    // look at each small regions (islands)
    for( int j = 0; j < nRegions; j++ )
      {
      if( j != largestId )
        {
        connect->InitializeSpecifiedRegionList();
        connect->AddSpecifiedRegion(j);
        connect->SetExtractionModeToSpecifiedRegions();

        island = connect->GetOutput();
#if (VTK_MAJOR_VERSION < 6)
        island->Update();
#endif
        island->BuildLinks();

        int ncells = island->GetNumberOfCells();
        // look at each point's LabelValue
        // if it has maximum #of not being "the Label",
        // set the non-label to all of the points on the island.
        // clean up neighborLabels first
        for( int ii = 0; ii < nLabels; ii++ )
          {
          neighborLabels[ii] = 0;
          }

        vtkIdType  npts;
        vtkIdType *pts;
        double     pt[3];
        // put borderIds into a list
        // Ids are from original surface
        vtkSmartPointer<vtkIdList> borderIds = vtkSmartPointer<vtkIdList>::New();
        for( int ii = 0; ii < ncells; ii++ )
          {
          island->GetCellPoints(ii, npts, pts);
          // decide if the point is on the border
          // by checking it's label
          for( int jj = 0; jj < npts; jj++ )
            {
            // get the label
            island->GetPoint(pts[jj], pt);
            int pid_orig = mesh->FindPoint(pt[0], pt[1], pt[2]);
            int label_jj = labelArray->GetTuple1(pid_orig);
            if( label_jj != label_i )
              {
              borderIds->InsertUniqueId(pid_orig);
              }
            }
          }
        // add neighbor labels
        for( int ii = 0; ii < borderIds->GetNumberOfIds(); ii++ )
          {
          int label_ii = labelArray->GetTuple1(borderIds->GetId(ii) );
          neighborLabels[label_ii - int(labelRange[0])] += 1;
          }

        // find the newLabel to be the one in neighborLabels with
        // maximum frequency
        int newLabel = label_i;
        for( int ii = 0; ii < nLabels; ii++ )
          {
          if( neighborLabels[ii] > neighborLabels[newLabel - int(labelRange[0])] )
            {
            newLabel = ii + labelRange[0];
            }
          }

        if( newLabel != label_i )
          {
          numChanges += 1;
          }
        // if found newLabel, assign it to all of the points of the island
        // on mesh
        for( int ii = 0; ii < ncells; ii++ )
          {
          island->GetCellPoints(ii, npts, pts);
          for( int jj = 0; jj < npts; jj++ )
            {
            // get each point
            island->GetPoint(pts[jj], pt);
            int pid_orig = mesh->FindPoint(pt[0], pt[1], pt[2]);
            labelArray->SetTuple1(pid_orig, newLabel);
            }
          }
        }
      }
    }
  delete [] neighborLabels;
  return numChanges;
}

// ----------------------------------------------------------------------------------------
// The function removes isolated points from the surface
int RemoveIsolatedPoints(vtkSmartPointer<vtkPolyData> mesh)
{
  // initialize number of changes made in this function
  int numChanges = 0;

  vtkDataArray *labelArray = mesh->GetPointData()->GetScalars();

  // mesh to get ready for GetPointCells
  mesh->BuildLinks();

  // label range
  double labelRange[2];
  labelArray->GetRange(labelRange);

  // number of labels
  int                           nLabels = int(labelRange[1] - labelRange[0] + 1);
  int *                         neighborLabels = new int[nLabels]; // keep neighbor labels
  vtkSmartPointer<vtkMaskLabel> mask = vtkSmartPointer<vtkMaskLabel>::New();
  for( int i = 0; i < nLabels; i++ )
    {
    int label_i = labelRange[0] + i;

    // mask out the label regions first
    // mask filter keeps the number of points
#if (VTK_MAJOR_VERSION < 6)
    mask->SetInput(mesh);
#else
    mask->SetInputData(mesh);
#endif
    mask->SetLabel(label_i);
    mask->Update();

    vtkPolyData *labelMask_i = mask->GetOutput();
    // labelMask_i->BuildLinks();
    int ncells = labelMask_i->GetNumberOfCells();
    // go throught the cells
    for( int j = 0; j < ncells; j++ )
      {
      vtkIdType  npts;
      vtkIdType *pts;
      double     pt[3];
      labelMask_i->GetCellPoints(j, npts, pts);
      // go through each point of the cell
      for( int jj = 0; jj < npts; jj++ )
        {
        // point Ids are different on labelMask_i
        // go by the physical location
        labelMask_i->GetPoint(pts[jj], pt);
        int pid_orig = mesh->FindPoint(pt[0], pt[1], pt[2]);
        // get the label of it through original Id
        int label_jj = labelArray->GetTuple1(pid_orig);

        // if center point has the same label as the mask label
        if( label_jj == label_i )
          {
          // go to its cells
          vtkSmartPointer<vtkIdList> cellIds = vtkSmartPointer<vtkIdList>::New();
          mesh->GetPointCells(pid_orig, cellIds);

          // check the labels of each point belongs to cellIds
          int nCells_jj = cellIds->GetNumberOfIds();

          vtkIdType *                pts_ii;
          vtkIdType                  npts_ii;
          vtkIdType                  cellId_ii;
          vtkSmartPointer<vtkIdList> neighborIds = vtkSmartPointer<vtkIdList>::New();
          // set up an array to keep neighbor labels
          for( int ii = 0; ii < nLabels; ii++ )
            {
            neighborLabels[ii] = 0;
            }
          for( int ii = 0; ii < nCells_jj; ii++ )
            {
            cellId_ii = cellIds->GetId(ii);
            mesh->GetCellPoints(cellId_ii, npts_ii, pts_ii);
            // go through neighbor cell's points
            // which are neighbor points of pid_orig
            for( int ni = 0; ni < npts_ii; ni++ )
              {
              // keep the neighbor Id if label is different with the center
              if( pts_ii[ni] != pid_orig )
                {
                neighborIds->InsertUniqueId(pts_ii[ni]);
                }
              }
            }
          // add neighbor labels
          // calculate number of neighbors having different labels with the center
          int nAliens = 0;
          for( int neighbor_i = 0; neighbor_i < neighborIds->GetNumberOfIds(); neighbor_i++ )
            {
            int label_n = labelArray->GetTuple1(neighborIds->GetId(neighbor_i) );
            if( label_n != label_i )
              {
              neighborLabels[label_n - int(labelRange[0])] += 1;
              nAliens += 1;
              }
            }

          // find the newLabel to be the one in neighborLabels with
          // maximum frequency
          int newLabel = label_i;
          for( int ii = 0; ii < nLabels; ii++ )
            {
            if( neighborLabels[ii] > neighborLabels[newLabel - int(labelRange[0])] )
              {
              newLabel = ii + int(labelRange[0]);
              }
            }

          int nNeighbors = neighborIds->GetNumberOfIds();

          // if surrounded by nCells_jj points that is different with mask label
          if( nNeighbors == nAliens )
            {
            labelArray->SetTuple1(pid_orig, newLabel);
            numChanges += 1;
            // std::cout<<"change happens at isolated points"<<std::endl;
            }
          }
        }
      }
    }
  delete [] neighborLabels;
  return numChanges;
}

// ----------------------------------------------------------------------------------------
// The function changes labels for sharp triangles on the border
// change the tip point's label
// there is not return for this function
void FlipSharpTriangles(vtkSmartPointer<vtkPolyData> mesh)
{
  mesh->BuildLinks();  // ready for GetPointCells

  vtkDataArray *labelArray = mesh->GetPointData()->GetScalars();

  // label range
  double labelRange[2];
  labelArray->GetRange(labelRange);

  // number of labels
  int  nLabels = int(labelRange[1] - labelRange[0] + 1);
  int *neighborLabels = new int[nLabels];   // keep neighbor labels
  int  iter = 0;
  int  change = 1;
  while( change > 0 && iter < 10 )
    {
    // clean up number of changes
    change = 0;
    // mask out each label with pure labels left
    vtkSmartPointer<vtkMaskLabel> mask = vtkSmartPointer<vtkMaskLabel>::New();
    for( int i = 0; i < nLabels; i++ )
      {
      int label_i = labelRange[0] + i;

      // mask out the label regions first
      // mask filter keeps the number of points
#if (VTK_MAJOR_VERSION < 6)
      mask->SetInput(mesh);
#else
      mask->SetInputData(mesh);
#endif
      mask->SetLabel(label_i);
      mask->LabelOnlyOn();
      mask->Update();

      vtkPolyData *mask_i = mask->GetOutput();
      mask_i->BuildLinks();
      // get point original Ids (on surface_in) for
      // tip points on borders
      int                        nCells = mask_i->GetNumberOfCells();
      vtkSmartPointer<vtkIdList> tipIds = vtkSmartPointer<vtkIdList>::New();
      // go through cells
      vtkIdType  npts;
      vtkIdType *pts;
      double     pt[3];
      for( int j = 0; j < nCells; j++ )
        {
        mask_i->GetCellPoints(j, npts, pts);
        // check neighbor cells of each edge
        vtkSmartPointer<vtkIdList> edgeCells0 = vtkSmartPointer<vtkIdList>::New();
        vtkSmartPointer<vtkIdList> edgeCells1 = vtkSmartPointer<vtkIdList>::New();
        vtkSmartPointer<vtkIdList> edgeCells2 = vtkSmartPointer<vtkIdList>::New();
        mask_i->GetCellEdgeNeighbors(j, pts[0], pts[1], edgeCells0);
        mask_i->GetCellEdgeNeighbors(j, pts[0], pts[2], edgeCells1);
        mask_i->GetCellEdgeNeighbors(j, pts[1], pts[2], edgeCells2);
        // keep sharp points original Id
        // the point with two edges having zero neighbor cells
        int pId_orig = -1;
        if( edgeCells0->GetNumberOfIds() == 0 && edgeCells1->GetNumberOfIds() == 0 )
          {
          // it is pts[0]
          mask_i->GetPoint(pts[0], pt);
          pId_orig = mesh->FindPoint(pt[0], pt[1], pt[2]);
          }
        if( edgeCells0->GetNumberOfIds() == 0 && edgeCells2->GetNumberOfIds() == 0 )
          {
          // it is pts[1]
          mask_i->GetPoint(pts[1], pt);
          pId_orig = mesh->FindPoint(pt[0], pt[1], pt[2]);
          }
        if( edgeCells1->GetNumberOfIds() == 0 && edgeCells2->GetNumberOfIds() == 0 )
          {
          // it is pts[2]
          mask_i->GetPoint(pts[2], pt);
          pId_orig = mesh->FindPoint(pt[0], pt[1], pt[2]);
          }
        // save the Id if found one
        if( pId_orig >= 0 )
          {
          tipIds->InsertUniqueId(pId_orig);
          }
        }

      // change labels for tip points
      // by finding the max number of neighbor labels
      int nTips = tipIds->GetNumberOfIds();
      for( int j = 0; j < nTips; j++ )
        {
        // clean up neighborLabels first
        for( int ii = 0; ii < nLabels; ii++ )
          {
          neighborLabels[ii] = 0;
          }
        int                        tId = tipIds->GetId(j);
        vtkSmartPointer<vtkIdList> cellIds = vtkSmartPointer<vtkIdList>::New();
        vtkSmartPointer<vtkIdList> neighborIds = vtkSmartPointer<vtkIdList>::New();
        mesh->GetPointCells(tId, cellIds);
        // go through neighbor points (cell points)
        // if its label is different with label_i
        // keep the label
        for( int jj = 0; jj < cellIds->GetNumberOfIds(); jj++ )
          {
          mesh->GetCellPoints(cellIds->GetId(jj), npts, pts);
          for( int cellp = 0; cellp < npts; cellp++ )
            {
            if( pts[cellp] != tId )
              {
              neighborIds->InsertUniqueId(pts[cellp]);
              }
            }
          }
        // add neighbor labels
        for( int neighbor_i = 0; neighbor_i < neighborIds->GetNumberOfIds(); neighbor_i++ )
          {
          int label_n = labelArray->GetTuple1(neighborIds->GetId(neighbor_i) );
          if( label_n != label_i )
            {
            neighborLabels[label_n - int(labelRange[0])] += 1;
            }
          }
        // find a new label for the tip point
        int newLabel = label_i;
        for( int ii = 0; ii < nLabels; ii++ )
          {
          if( neighborLabels[ii] > neighborLabels[newLabel - int(labelRange[0])] )
            {
            newLabel = ii + labelRange[0];
            }
          }

        if( newLabel != label_i )
          {
          // set new label to the tip point
          labelArray->SetTuple1(tId, newLabel);
          change += 1;
          }
        }
      }
    iter += 1;
    }

  delete [] neighborLabels;
}
