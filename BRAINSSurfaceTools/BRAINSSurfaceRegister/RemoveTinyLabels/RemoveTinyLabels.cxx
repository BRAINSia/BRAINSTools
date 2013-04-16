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
#include "vtkDataArray.h"
#include "vtkPolyDataConnectivityFilter.h"

#include "vtkMaskLabel.h"
#include "RemoveTinyLabelsCLP.h"
#include "BRAINSvtkV6Compat.h"

int main( int argc, char * argv[] )
{
  PARSE_ARGS;

  std::cout << "-----------------------------------------------" << std::endl;
  std::cout << "Input Surface: " << std::endl;
  std::cout << inputSurfaceFile << std::endl;
  std::cout << "Output Surface: " << std::endl;
  std::cout << outputSurfaceFile << std::endl;
  std::cout << "Remove Labels: " << std::endl;
  for( unsigned int i = 0; i < labelList.size(); i++ )
    {
    std::cout << labelList[i] << "; ";
    }
  std::cout << std::endl;
  std::cout << "-----------------------------------------------" << std::endl;

  // Create all of the classes we will need
  vtkSmartPointer<vtkPolyDataReader> reader =
    vtkSmartPointer<vtkPolyDataReader>::New();
  vtkSmartPointer<vtkPolyData> island =
    vtkSmartPointer<vtkPolyData>::New();
  vtkSmartPointer<vtkMaskLabel> mask =
    vtkSmartPointer<vtkMaskLabel>::New();
  vtkSmartPointer<vtkPolyDataConnectivityFilter> connect =
    vtkSmartPointer<vtkPolyDataConnectivityFilter>::New();
  vtkSmartPointer<vtkPolyDataWriter> writer =
    vtkSmartPointer<vtkPolyDataWriter>::New();

  // Define all of the variables
  unsigned int nAbsentLabels = labelList.size();
  // unsigned int *absentLabel = new unsigned int[50];
  unsigned int *neighborLabels = new unsigned int[50];
  int           ncells;
  vtkIdType *   pts = 0;
  vtkIdType     npts = 0;
  vtkIdType     pid_orig;
  double *      pt;
  unsigned int  newLabel;
  unsigned int  label, label_jj;

  // read the label surface
  reader->SetFileName(inputSurfaceFile.c_str() );
  reader->Update();
  vtkSmartPointer<vtkPolyData> surface_in = reader->GetOutput();

  vtkDataArray *labelArray = surface_in->GetPointData()->GetScalars();

  std::string labelName = labelArray->GetName();

  if( (labelArray == NULL) || (labelName != "LabelValue") )
    {
    std::cerr << "There is no labelarray on the input surface. ";
    std::cerr << "Quit." << std::endl;
    return 1;
    }
  // go through each label in absentLabel
  for( unsigned int i = 0; i < nAbsentLabels; i++ )
    {
    label = labelList[i];

    // std::cout<<"remove label: "<<label<<std::endl;

    // analyze each label in the list
    BRAINSvtkV6_SetInputData( mask, reader->GetOutput() );
    mask->SetLabel(label);
    mask->Update();

    // replace the labels on mask->output if it is not null
    if( mask->GetOutput()->GetNumberOfCells() )
      {
      BRAINSvtkV6_SetInputData( connect, mask->GetOutput() );
      connect->SetExtractionModeToAllRegions();
      connect->Update();

      int nRegions = connect->GetNumberOfExtractedRegions();
      // look at each region
      for( int j = 0; j < nRegions; j++ )
        {
        connect->InitializeSpecifiedRegionList();
        connect->AddSpecifiedRegion(j);
        connect->SetExtractionModeToSpecifiedRegions();

        connect->Update();
        island = connect->GetOutput();
        island->BuildLinks();

        ncells = island->GetNumberOfCells();
        // std::cout<<"label "<<label<<" has "<<ncells<<" cells."<<std::endl;
        // look at each point's LabelValue
        // find out another label with maximum #of points having "not the Label",
        // set the non-label to all of the points on the island.
        // clean up neighborLabels first
        // assume the input surface has no more than 50 labels
        for( unsigned int ii = 0; ii < 50; ii++ )
          {
          neighborLabels[ii] = 0;
          }
        for( int ii = 0; ii < ncells; ii++ )
          {
          island->GetCellPoints(ii, npts, pts);
          for( int jj = 0; jj < npts; jj++ )
            {
            // get each point
            pt = island->GetPoint(pts[jj]);
            pid_orig = surface_in->FindPoint(pt[0], pt[1], pt[2]);
            label_jj = labelArray->GetTuple1(pid_orig);
            if( label_jj != label )
              {
              neighborLabels[label_jj] += 1;
              }
            }
          }

        // find the newLabel to be the one in neighborLabels with
        // maximum frequency
        newLabel = label;
        for( unsigned int ii = 0; ii < 50; ii++ )
          {
          if( neighborLabels[ii] > neighborLabels[newLabel] )
            {
            newLabel = ii;
            }
          }
        // assign newLabel to all of the points of the island
        // on surface_in
        for( int ii = 0; ii < ncells; ii++ )
          {
          island->GetCellPoints(ii, npts, pts);
          for( int jj = 0; jj < npts; jj++ )
            {
            // get each point
            pt = island->GetPoint(pts[jj]);
            pid_orig = surface_in->FindPoint(pt[0], pt[1], pt[2]);
            labelArray->SetTuple1(pid_orig, newLabel);
            }
          }
        }
      }
    }

  BRAINSvtkV6_SetInputData(writer, surface_in);
  writer->SetFileName(outputSurfaceFile.c_str() );
  writer->Update();

  delete [] neighborLabels;
  return 0;
}
