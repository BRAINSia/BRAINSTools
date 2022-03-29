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

  Program:   Visualization Toolkit
  Module:    $RCSfile: vtkStructuredGrid.cxx,v $

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
#include "vtkObjectFactory.h"
#include "vtkInformation.h"
#include "vtkStructuredGrid.h"
#include "vtkUnstructuredGrid.h"
#include "vtkTalairachGrid.h"
#include "vtkStructuredGridWriter.h"
#include "itkMacro.h" //Needed for nullptr
#define PR(x)                                                                                                          \
  std::cout << #x " = " << x << "\n"; // a simple print macro for
                                      // use when debugging

vtkStandardNewMacro(vtkTalairachGrid);

void
vtkTalairachGrid::Initialize()
{
  boundingBoxGrid->Initialize();
  boundingBoxGrid->SetDimensions(0, 0, 0);
  talairachGrid->Initialize();
  talairachGrid->SetDimensions(0, 0, 0);
}

int *
vtkTalairachGrid::GetBoundingBoxDimensions()
{
  return (this->boundingBoxGrid)->GetDimensions();
}

int *
vtkTalairachGrid::GetTalairachGridDimensions()
{
  return (this->talairachGrid)->GetDimensions();
}

void
vtkTalairachGrid::SetACPoint(double ac[3])
{
  ACPoint[0] = ac[0];
  ACPoint[1] = ac[1];
  ACPoint[2] = ac[2];
}

void
vtkTalairachGrid::SetPCPoint(double pc[3])
{
  PCPoint[0] = pc[0];
  PCPoint[1] = pc[1];
  PCPoint[2] = pc[2];
}

void
vtkTalairachGrid::SetIRPPoint(double irp[3])
{
  IRPPoint[0] = irp[0];
  IRPPoint[1] = irp[1];
  IRPPoint[2] = irp[2];
}

void
vtkTalairachGrid::SetSLAPoint(double sla[3])
{
  SLAPoint[0] = sla[0];
  SLAPoint[1] = sla[1];
  SLAPoint[2] = sla[2];
}

double *
vtkTalairachGrid::GetACPoint()
{
  return this->ACPoint;
}

double *
vtkTalairachGrid::GetPCPoint()
{
  return this->PCPoint;
}

double *
vtkTalairachGrid::GetIRPPoint()
{
  return this->IRPPoint;
}

double *
vtkTalairachGrid::GetSLAPoint()
{
  return this->SLAPoint;
}

vtkStructuredGrid *
vtkTalairachGrid::GetTalairachGrid()
{
  std::cout << "grid pnts in talgrid::GetTalairachGrid:: " << this->talairachGrid->GetNumberOfPoints() << std::endl;

  return this->talairachGrid;
}

vtkStructuredGrid *
vtkTalairachGrid::GetBoundingBoxGrid()
{
  return this->boundingBoxGrid;
}

void
vtkTalairachGrid::SetTalairachGrid(vtkStructuredGrid * grid)
{
  // this->talairachGrid = grid;
  // talairachGrid->DeepCopy(grid);
  // talairachGridPoints->DeepCopy(grid->GetPoints());
  std::cout << "grid pnts in talgrid.cxx1: " << grid->GetNumberOfPoints() << std::endl;

  // this->talairachGrid = grid;
  this->talairachGrid->DeepCopy(grid);
  std::cout << "grid pnts in talgrid.cxx2: " << this->talairachGrid->GetNumberOfPoints() << std::endl;
}

void
vtkTalairachGrid::SetBoundingBoxGrid(vtkStructuredGrid * box)
{
  // this->boundingBoxGrid = box;
  boundingBoxGrid->DeepCopy(box);
}

vtkPoints *
vtkTalairachGrid::GetTalairachGridPoints()
{
  return this->talairachGridPoints;
}

vtkPoints *
vtkTalairachGrid::GetBoundingBoxGridPoints()
{
  return this->boundingBoxGridPoints;
}

void
vtkTalairachGrid::EstablishBoundingBoxGrid()
{
  vtkPoints * points = vtkPoints::New();
  int *       dim = GetBoundingBoxDimensions();

  points->Allocate(dim[0] * dim[1] * dim[2]);

  double              pnt[3];
  std::vector<double> talPnt;

  /* First Row */
  pnt[0] = IRPPoint[0];
  pnt[1] = IRPPoint[1];
  pnt[2] = IRPPoint[2];
  boundingBoxGridPoints->InsertPoint(0, pnt);

  pnt[0] = ACPoint[0];
  pnt[1] = IRPPoint[1];
  pnt[2] = IRPPoint[2];
  boundingBoxGridPoints->InsertPoint(1, pnt);

  pnt[0] = SLAPoint[0];
  pnt[1] = IRPPoint[1];
  pnt[2] = IRPPoint[2];
  boundingBoxGridPoints->InsertPoint(2, pnt);

  /* Second Row */
  pnt[0] = IRPPoint[0];
  pnt[1] = ACPoint[1];
  pnt[2] = IRPPoint[2];
  boundingBoxGridPoints->InsertPoint(3, pnt);

  pnt[0] = ACPoint[0];
  pnt[1] = ACPoint[1];
  pnt[2] = IRPPoint[2];
  boundingBoxGridPoints->InsertPoint(4, pnt);

  pnt[0] = SLAPoint[0];
  pnt[1] = ACPoint[1];
  pnt[2] = IRPPoint[2];
  boundingBoxGridPoints->InsertPoint(5, pnt);

  /* Third Row */
  pnt[0] = IRPPoint[0];
  pnt[1] = PCPoint[1];
  pnt[2] = IRPPoint[2];
  boundingBoxGridPoints->InsertPoint(6, pnt);

  pnt[0] = ACPoint[0];
  pnt[1] = PCPoint[1];
  pnt[2] = IRPPoint[2];
  boundingBoxGridPoints->InsertPoint(7, pnt);

  pnt[0] = SLAPoint[0];
  pnt[1] = PCPoint[1];
  pnt[2] = IRPPoint[2];
  boundingBoxGridPoints->InsertPoint(8, pnt);

  /* Fourth Row */
  pnt[0] = IRPPoint[0];
  pnt[1] = SLAPoint[1];
  pnt[2] = IRPPoint[2];
  boundingBoxGridPoints->InsertPoint(9, pnt);

  pnt[0] = ACPoint[0];
  pnt[1] = SLAPoint[1];
  pnt[2] = IRPPoint[2];
  boundingBoxGridPoints->InsertPoint(10, pnt);

  pnt[0] = SLAPoint[0];
  pnt[1] = SLAPoint[1];
  pnt[2] = (IRPPoint[2]);
  boundingBoxGridPoints->InsertPoint(11, pnt);
  /**** End First Slice ***/

  /***Start Second Slice ***/
  /* First Slice */
  pnt[0] = IRPPoint[0];
  pnt[1] = IRPPoint[1];
  pnt[2] = ACPoint[2];
  boundingBoxGridPoints->InsertPoint(12, pnt);

  pnt[0] = ACPoint[0];
  pnt[1] = IRPPoint[1];
  pnt[2] = ACPoint[2];
  boundingBoxGridPoints->InsertPoint(13, pnt);

  pnt[0] = SLAPoint[0];
  pnt[1] = IRPPoint[1];
  pnt[2] = ACPoint[2];
  boundingBoxGridPoints->InsertPoint(14, pnt);

  /* Second Row */
  pnt[0] = IRPPoint[0];
  pnt[1] = ACPoint[1];
  pnt[2] = ACPoint[2];
  boundingBoxGridPoints->InsertPoint(15, pnt);

  pnt[0] = ACPoint[0];
  pnt[1] = ACPoint[1];
  pnt[2] = ACPoint[2];
  boundingBoxGridPoints->InsertPoint(16, pnt);

  pnt[0] = SLAPoint[0];
  pnt[1] = ACPoint[1];
  pnt[2] = ACPoint[2];
  boundingBoxGridPoints->InsertPoint(17, pnt);

  /* Third Row */
  pnt[0] = IRPPoint[0];
  pnt[1] = PCPoint[1];
  pnt[2] = ACPoint[2];
  boundingBoxGridPoints->InsertPoint(18, pnt);

  pnt[0] = ACPoint[0];
  pnt[1] = PCPoint[1];
  pnt[2] = ACPoint[2];
  boundingBoxGridPoints->InsertPoint(19, pnt);

  pnt[0] = SLAPoint[0];
  pnt[1] = PCPoint[1];
  pnt[2] = ACPoint[2];
  boundingBoxGridPoints->InsertPoint(20, pnt);

  /* Fourth Row */
  pnt[0] = IRPPoint[0];
  pnt[1] = SLAPoint[1];
  pnt[2] = ACPoint[2];
  boundingBoxGridPoints->InsertPoint(21, pnt);

  pnt[0] = ACPoint[0];
  pnt[1] = SLAPoint[1];
  pnt[2] = ACPoint[2];
  boundingBoxGridPoints->InsertPoint(22, pnt);

  pnt[0] = SLAPoint[0];
  pnt[1] = SLAPoint[1];
  pnt[2] = ACPoint[2];
  boundingBoxGridPoints->InsertPoint(23, pnt);
  /***End of Second Slice ***/

  /*** Third Slice ***/
  /* First Row */
  pnt[0] = IRPPoint[0];
  pnt[1] = IRPPoint[1];
  pnt[2] = (SLAPoint[2]);
  boundingBoxGridPoints->InsertPoint(24, pnt);

  pnt[0] = ACPoint[0];
  pnt[1] = IRPPoint[1];
  pnt[2] = SLAPoint[2];
  boundingBoxGridPoints->InsertPoint(25, pnt);

  pnt[0] = SLAPoint[0];
  pnt[1] = IRPPoint[1];
  pnt[2] = SLAPoint[2];
  boundingBoxGridPoints->InsertPoint(26, pnt);

  /* Second Row */
  pnt[0] = IRPPoint[0];
  pnt[1] = ACPoint[1];
  pnt[2] = SLAPoint[2];
  boundingBoxGridPoints->InsertPoint(27, pnt);

  pnt[0] = ACPoint[0];
  pnt[1] = ACPoint[1];
  pnt[2] = SLAPoint[2];
  boundingBoxGridPoints->InsertPoint(28, pnt);

  pnt[0] = SLAPoint[0];
  pnt[1] = ACPoint[1];
  pnt[2] = SLAPoint[2];
  boundingBoxGridPoints->InsertPoint(29, pnt);

  /* Third Row */
  pnt[0] = IRPPoint[0];
  pnt[1] = PCPoint[1];
  pnt[2] = SLAPoint[2];
  boundingBoxGridPoints->InsertPoint(30, pnt);

  pnt[0] = ACPoint[0];
  pnt[1] = PCPoint[1];
  pnt[2] = SLAPoint[2];
  boundingBoxGridPoints->InsertPoint(31, pnt);

  pnt[0] = SLAPoint[0];
  pnt[1] = PCPoint[1];
  pnt[2] = SLAPoint[2];
  boundingBoxGridPoints->InsertPoint(32, pnt);

  /* Fourth Row */
  pnt[0] = IRPPoint[0];
  pnt[1] = SLAPoint[1];
  pnt[2] = SLAPoint[2];
  boundingBoxGridPoints->InsertPoint(33, pnt);

  pnt[0] = ACPoint[0];
  pnt[1] = SLAPoint[1];
  pnt[2] = SLAPoint[2];
  boundingBoxGridPoints->InsertPoint(34, pnt);

  pnt[0] = SLAPoint[0];
  pnt[1] = SLAPoint[1];
  pnt[2] = SLAPoint[2];
  boundingBoxGridPoints->InsertPoint(35, pnt);
}

void
vtkTalairachGrid::EstablishTalairachGrid()
{
  // int* dim = GetTalairachGridDimensions();
  // vtkPoints *points = vtkPoints::New();
  // points->Allocate(dim[0]*dim[1]*dim[2]);

  // short x_alloc = static_cast<short>(dim[0]);
  // short y_alloc = static_cast<short>(dim[1]);
  // short z_alloc = static_cast<short>(dim[2]);

  // points->Allocate(x_alloc * y_alloc * z_alloc);

  int offset = 0;

  int i = 0;
  int j = 0;
  int k = 0;

  double pnt[3] = { 0.0, 0.0, 0.0 };

  double iVal = 0.0;
  double jVal = 0.0;
  double kVal = 0.0;

  double x0 = (IRPPoint[0]);
  double x4 = (PCPoint[0] + 1.0);
  double x8 = (SLAPoint[0]);

  // double z0 = floor(SLAPoint[2] + 1.0);
  double z0 = (IRPPoint[2]);

  double z8 = (PCPoint[2]);
  // double z12 = floor(IRPPoint[2] + 1.0);
  double z12 = (SLAPoint[2]);

  // double z0 = floor((z8 + (zIrp - z8) / 4.0 * 6.0));
  for (k = 0; k < 15; k++)
  {
    switch (k)
    {
        /*
              case 0: kVal = z0;
                      break;
              case 1: kVal = floor((z0 + (z8 - z0) / 8.0));
                      break;
              case 2: kVal = floor((z0 + (z8 - z0) / 8.0 * 2.0));
                      break;
              case 3: kVal = floor((z0 + (z8 - z0) / 8.0 * 3.0));
                      break;
              case 4: kVal = floor((z0 + (z8 - z0) / 8.0 * 4.0));
                      break;
              case 5: kVal = floor((z0 + (z8 - z0) / 8.0 * 5.0));
                      break;
              case 6: kVal = floor((z0 + (z8 - z0) / 8.0 * 6.0));
                      break;
              case 7: kVal = floor((z0 + (z8 - z0) / 8.0 * 7.0));
                      break;
              case 8: kVal = z8;
                      break;
              case 9: kVal = floor((z8 + (z12 - z8) / 4.0));
                      break;
              case 10: kVal = floor((z8 + (z12 - z8) / 4.0 * 2.0));
                       break;
              case 11: kVal = floor((z8 + (z12 - z8) / 4.0 * 3.0));
                       break;
              case 12: kVal = z12;
                       break;
              case 13: kVal = floor((z8 + (z12 - z8) / 4.0 * 5.0));
                       break;
              case 14: kVal = floor((z8 + (z12 - z8) / 4.0 * 6.0));
                       break;
              default: kVal = z0;
                       break;
        */

      case 0:
      {
        kVal = ((z8 + (z0 - z8) / 4.0 * 6.0));
      }
      break;
      case 1:
      {
        kVal = ((z8 + (z0 - z8) / 4.0 * 5.0));
      }
      break;
      case 2:
      {
        kVal = ((z8 + (z0 - z8) / 4.0 * 4.0));
      }
      break;
      case 3:
      {
        kVal = ((z8 + (z0 - z8) / 4.0 * 3.0));
      }
      break;
      case 4:
      {
        kVal = ((z8 + (z0 - z8) / 4.0 * 2.0));
      }
      break;
      case 5:
      {
        kVal = ((z8 + (z0 - z8) / 4.0));
      }
      break;
      case 6:
      {
        kVal = z8;
      }
      break;
      case 7:
      {
        kVal = ((z8 + (z12 - z8) / 8.0));
      }
      break;
      case 8:
      {
        kVal = ((z8 + (z12 - z8) / 8.0 * 2.0));
      }
      break;
      case 9:
      {
        kVal = ((z8 + (z12 - z8) / 8.0 * 3.0));
      }
      break;
      case 10:
      {
        kVal = ((z8 + (z12 - z8) / 8.0 * 4.0));
      }
      break;
      case 11:
      {
        kVal = ((z8 + (z12 - z8) / 8.0 * 5.0));
      }
      break;
      case 12:
      {
        kVal = ((z8 + (z12 - z8) / 8.0 * 6.0));
      }
      break;
      case 13:
      {
        kVal = ((z8 + (z12 - z8) / 8.0 * 7.0));
      }
      break;
      case 14:
      {
        kVal = z12;
      }
      break;
      default:
      {
        kVal = z0;
      }
      break;
    }

    jVal = 0.0;
    for (j = 0; j < 12; j++)
    {
      switch (j)
      {
        case 0:
        {
          jVal = IRPPoint[1];
        }
        break;
        case 1:
        {
          jVal = (IRPPoint[1] + (ACPoint[1] - IRPPoint[1]) / 4.0);
        }
        break;
        case 2:
        {
          jVal = ((IRPPoint[1] + (ACPoint[1] - IRPPoint[1]) / 4.0 * 2.0));
        }
        break;
        case 3:
        {
          jVal = ((IRPPoint[1] + (ACPoint[1] - IRPPoint[1]) / 4.0 * 3.0));
        }
        break;
        case 4:
        {
          jVal = ACPoint[1];
        }
        break;
        case 5:
        {
          jVal = (PCPoint[1] + ((ACPoint[1] - PCPoint[1]) / 3.0) * 2.0);
        }
        break;
        case 6:
        {
          jVal = (PCPoint[1] + ((ACPoint[1] - PCPoint[1]) / 3.0));
        }
        break;
        case 7:
        {
          jVal = PCPoint[1];
        }
        break;
        case 8:
        {
          jVal = ((PCPoint[1] + (SLAPoint[1] - PCPoint[1]) / 4.0));
        }
        break;
        case 9:
        {
          jVal = ((PCPoint[1] + (SLAPoint[1] - PCPoint[1]) / 4.0 * 2.0));
        }
        break;
        case 10:
        {
          jVal = ((PCPoint[1] + (SLAPoint[1] - PCPoint[1]) / 4.0 * 3.0));
        }
        break;
        case 11:
        {
          jVal = SLAPoint[1];
        }
        break;
        default:
        {
          jVal = IRPPoint[1];
        }
        break;
      }

      iVal = 0.0;
      for (i = 0; i < 9; i++)
      {
        switch (i)
        {
          case 0:
          {
            iVal = x0;
          }
          break;
          case 1:
          {
            iVal = ((x4 + (x0 - x4) / 4.0 * 3.0));
          }
          break;
          case 2:
          {
            iVal = ((x4 + (x0 - x4) / 4.0 * 2.0));
          }
          break;
          case 3:
          {
            iVal = ((x4 + (x0 - x4) / 4.0));
          }
          break;
          case 4:
          {
            iVal = x4;
          }
          break;
          case 5:
          {
            iVal = ((x4 + (x8 - x4) / 4.0));
          }
          break;
          case 6:
          {
            iVal = ((x4 + (x8 - x4) / 4.0 * 2.0));
          }
          break;
          case 7:
          {
            iVal = ((x4 + (x8 - x4) / 4.0 * 3.0));
          }
          break;
          case 8:
          {
            iVal = x8;
          }
          break;
          default:
          {
            iVal = x0;
          }
          break;
        }

        pnt[0] = iVal;
        pnt[1] = jVal;
        pnt[2] = kVal;

        // std::vector<double> talPnt;
        // talPnt = ConvertPixelPointToTalairachPoint(pnt);
        // pnt[0] = talPnt[0];
        // pnt[1] = talPnt[1];
        // pnt[2] = talPnt[2];

        talairachGridPoints->InsertPoint(offset, pnt);
        offset++;
      }
    }
  }
}

void
vtkTalairachGrid::Update()
{
  boundingBoxGrid->SetPoints(boundingBoxGridPoints);
  talairachGrid->SetPoints(talairachGridPoints);
}

void
vtkTalairachGrid::WriteTalairachGrid(std::string filename)
{
  vtkStructuredGridWriter * gridWriter = vtkStructuredGridWriter::New();

  gridWriter->SetFileName(filename.c_str());
#if (VTK_MAJOR_VERSION < 6)
  gridWriter->SetInput(talairachGrid);
#else
  gridWriter->SetInputData(talairachGrid);
#endif
  gridWriter->Write();
}

void
vtkTalairachGrid::WriteBoundingBoxGrid(std::string filename)
{
  vtkStructuredGridWriter * gridWriter = vtkStructuredGridWriter::New();

  gridWriter->SetFileName(filename.c_str());
#if (VTK_MAJOR_VERSION < 6)
  gridWriter->SetInput(boundingBoxGrid);
#else
  gridWriter->SetInputData(boundingBoxGrid);
#endif

  gridWriter->Write();
}

void
vtkTalairachGrid::PrintSelf(ostream & os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);

  int boxDim[3];
  int gridDim[3];
  this->GetBoundingBoxGrid()->GetDimensions(boxDim);
  this->GetTalairachGrid()->GetDimensions(gridDim);

  os << indent << "Box Dimensions: (" << boxDim[0] << ", " << boxDim[1] << ", " << boxDim[2] << ")\n";

  os << indent << "Grid Dimensions: (" << gridDim[0] << ", " << gridDim[1] << ", " << gridDim[2] << ")\n";

  os << indent << "AC Point: " << this->GetACPoint() << "\n";
  os << indent << "PC Point: " << this->GetPCPoint() << "\n";
  os << indent << "SLA Point: " << this->GetSLAPoint() << "\n";
  os << indent << "IRP Point: " << this->GetIRPPoint() << "\n";

  const int * extent = this->Extent;
  os << indent << "Extent: " << extent[0] << ", " << extent[1] << ", " << extent[2] << ", " << extent[3] << ", "
     << extent[4] << ", " << extent[5] << endl;

  os << ")\n";
}

vtkTalairachGrid::vtkTalairachGrid()
{
  /* Initialize the 4 defining points to {0, 0, 0} at start */
  ACPoint[0] = 0;
  ACPoint[1] = 0;
  ACPoint[2] = 0;

  PCPoint[0] = 0;
  PCPoint[1] = 0;
  PCPoint[2] = 0;

  SLAPoint[0] = 0;
  SLAPoint[1] = 0;
  SLAPoint[2] = 0;

  IRPPoint[0] = 0;
  IRPPoint[1] = 0;
  IRPPoint[2] = 0;

  /* Instantiate the two grid representations */
  boundingBoxGridPoints = vtkPoints::New();
  talairachGridPoints = vtkPoints::New();

  /* Set default grid dimensions */
  boundingBoxGrid = vtkStructuredGrid::New();
  int dimensions[3];
  dimensions[0] = 3;
  dimensions[1] = 4;
  dimensions[2] = 3;
  (this->boundingBoxGrid)->SetDimensions(dimensions);

  talairachGrid = vtkStructuredGrid::New();
  dimensions[0] = 9;
  dimensions[1] = 12;
  dimensions[2] = 15;
  (this->talairachGrid)->SetDimensions(dimensions);

  /* Set the extent and information for the grid object */
  int extent[6] = { 0, -1, 0, -1, 0, -1 };
  memcpy((this->Extent), extent, 6 * sizeof(int));

  (this->Information)->Set(vtkDataObject::DATA_EXTENT_TYPE(), VTK_3D_EXTENT);
  (this->Information)->Set(vtkDataObject::DATA_EXTENT(), (this->Extent), 6);
}

vtkTalairachGrid::~vtkTalairachGrid()
{
  if ((this->boundingBoxGrid))
  {
    (this->boundingBoxGrid)->Delete();
    (this->boundingBoxGrid) = nullptr;
  }
  if ((this->talairachGrid))
  {
    (this->talairachGrid)->Delete();
    (this->talairachGrid) = nullptr;
  }
  if ((this->boundingBoxGridPoints))
  {
    (this->boundingBoxGridPoints)->Delete();
    (this->boundingBoxGridPoints) = nullptr;
  }
  if ((this->talairachGridPoints))
  {
    (this->talairachGridPoints)->Delete();
    (this->talairachGridPoints) = nullptr;
  }
}

std::vector<double>
vtkTalairachGrid::ConvertTalairachPointToPixelPoint(double * talPoint)
{
  /* voxel point that is to be generated from talairach point */
  std::vector<double> voxelPoint;

  voxelPoint.resize(3);

  /* double is assigned the y coordinate of the point */
  double distance = talPoint[1];
  if (distance >= 0.0)
  {
    voxelPoint[1] = (distance / 69.0) * (SLAPoint[1] - ACPoint[1]) + ACPoint[1];
  }
  else if (distance >= -24.0)
  {
    voxelPoint[1] = (distance / 24.0) * (ACPoint[1] - PCPoint[1]) + ACPoint[1];
  }
  else
  {
    voxelPoint[1] = (distance + 24.0) / (102.0 - 24.0) * (PCPoint[1] - IRPPoint[1]) + PCPoint[1];
  }

  /* double is assigned the x coordinates of the point */
  distance = talPoint[0];
  if (distance >= 0.0)
  {
    voxelPoint[0] = (PCPoint[0] + ((distance / 68.0) * fabs((SLAPoint[0] - PCPoint[0]))));
  }
  else
  {
    voxelPoint[0] = (PCPoint[0] + ((distance / 68.0) * fabs((IRPPoint[0] - PCPoint[0]))));
  }

  /* double is assigned the z coordinates of the point */
  distance = talPoint[2];
  if (distance >= 0.0)
  {
    voxelPoint[2] = ((distance / 75.0) * (SLAPoint[2] - PCPoint[2]) + PCPoint[2]);
  }
  else
  {
    voxelPoint[2] = ((distance / 43.0) * (PCPoint[2] - IRPPoint[2]) + PCPoint[2]);
  }

  return voxelPoint;
}

std::vector<double>
vtkTalairachGrid::ConvertPixelPointToTalairachPoint(double * voxelPoint)
{
  /* talairach point that is to be generated from voxel point */
  std::vector<double> talPoint;

  talPoint.resize(3);

  /* Y axis */
  if (voxelPoint[1] <= ACPoint[1])
  {
    std::cout << "Y - Greater than AC" << std::endl;
    talPoint[1] = (voxelPoint[1] - ACPoint[1]) / static_cast<float>(SLAPoint[1] - ACPoint[1]) * (102.0 - 24.0) - 24.0;
  }
  else if (voxelPoint[1] <= PCPoint[1])
  {
    std::cout << "Y - Greater than PC" << std::endl;
    talPoint[1] = (voxelPoint[1] - ACPoint[1]) / static_cast<float>(ACPoint[1] - PCPoint[1]) * 24.0;
  }
  else
  {
    std::cout << "Y - Less than AC and PC" << std::endl;
    talPoint[1] = (voxelPoint[1] - PCPoint[1]) / static_cast<float>(PCPoint[1] - IRPPoint[1]) * 69.0;
  }

  /* X axis */
  if (voxelPoint[0] >= PCPoint[0])
  {
    std::cout << "X - Greater than PC" << std::endl;
    talPoint[0] = (voxelPoint[0] - PCPoint[0]) / fabs((SLAPoint[0] - PCPoint[0])) * 68.0;
  }
  else
  {
    std::cout << "X - Less than PC" << std::endl;
    talPoint[0] = (voxelPoint[0] - PCPoint[0]) / fabs((IRPPoint[0] - PCPoint[0])) * 68.0;
  }

  /* Z axis */
  if (voxelPoint[2] > PCPoint[2])
  {
    std::cout << "Z - Greater than PC" << std::endl;
    talPoint[2] = (voxelPoint[2] - PCPoint[2]) / static_cast<float>(SLAPoint[2] - PCPoint[2]) * 75.0;
  }
  else
  {
    std::cout << "Z - Less than PC" << std::endl;
    talPoint[2] = (voxelPoint[2] - PCPoint[2]) / static_cast<float>(PCPoint[2] - IRPPoint[2]) * 43.0;
  }

  std::cout << "Tal Point : " << talPoint[0] << " " << talPoint[1] << " " << talPoint[2] << std::endl;

  return talPoint;
}
