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

#include "itkDtiTrackingFilterBase.h"
#include "algo.h"
#include "GtractTypes.h"
#include "gtractCommonWin32.h"

#include <map>
#include <string>

VTKFiberListType ConvertFiberToVTK(FiberListType & fiberList, int order, AnisotropyImageType::Pointer anisotropyImage)
{
  VTKFiberListType vtkFiberList;

  FiberListType::iterator fiberIt;

  for( fiberIt = fiberList.begin(); fiberIt != fiberList.end(); fiberIt++ )
    {
    FiberType      fiber = *fiberIt;
    vtkPoints *    points = vtkPoints::New();
    vtkFloatArray *scalars = vtkFloatArray::New();
    scalars->SetName("Anisotropy");
    int N = 0;
    if( order == 1 )
      {
      fiber.reverse();
      }

    while( !fiber.empty() )
      {
      FiberPointType       fp = fiber.front();
      itk::Point<float, 3> p;
      anisotropyImage->TransformContinuousIndexToPhysicalPoint(fp.m_Point, p);
      points->InsertPoint(N, p[0], p[1], p[2]);
      scalars->InsertTuple1(N, 1 - fp.m_AI);
      fiber.pop_front();
      N++;
      }

    vtkCellArray *line = vtkCellArray::New();
    line->InsertNextCell(N);
    for( int i = 0; i < N; i++ )
      {
      line->InsertCellPoint(i);
      }

    vtkPolyData *data = vtkPolyData::New();
    data->SetPoints(points);
    data->SetLines(line);
    data->GetPointData()->SetScalars(scalars);
    vtkFiberList.push_back(data);
    }

  return vtkFiberList;
}

VTKFiberListType MergeFibers(VTKFiberListType & fiberList1, VTKFiberListType & fiberList2)
{
  VTKFiberListType mergedFiberList;

  if( fiberList1.size() == 0 )
    {
    mergedFiberList = fiberList2;
    }
  else if( fiberList2.size() == 0 )
    {
    mergedFiberList = fiberList1;
    }
  else
    {
    std::cout << std::endl;
    std::cout << "Merging two fiber groups...." << std::endl;
    clock_t start, finish;
    start = clock();
    // ////////////////////////////////////////////////////////////////////////
    // merge fiber
    // ////////////////////////////////////////////////////////////////////////
    int     N1 = fiberList1.size();
    int     N2 = fiberList2.size();
    TVector dis1(N1);
    TVector dis2(N2);
    MinimumDistanceBetweenFiberGroups(fiberList1, fiberList2, dis1);
    MinimumDistanceBetweenFiberGroups(fiberList2, fiberList1, dis2);
    int disC, merC, index;
    disC = merC = 0;
    VTKFiberListType tempFiberList;
    tempFiberList.clear();
    double lowerT = 0.1;
    double upperT = 4;
    double mean, sd;
    double sum = 0;
    for( int i = 0; i < N1; i++ )
      {
      sum += dis1[i];
      }
    mean = sum / N1;
    sum = 0;
    for( int i = 0; i < N1; i++ )
      {
      sum += vcl_pow(dis1[i] - mean, 2.0);
      }
    sd = vcl_sqrt( sum / ( N1 - 1 ) );
    upperT = mean + 3 * sd;
    std::cout << "Mis-match threshold1: " << upperT << " Mean " << mean << std::endl;
    index = 0;
    while( !fiberList1.empty() )
      {
      if( dis1[index] < upperT )
        {
        mergedFiberList.push_back( fiberList1.front() );
        }
      else
        {
        disC++;
        }
      fiberList1.pop_front();
      index++;
      }

    sum = 0;
    for( int i = 0; i < N2; i++ )
      {
      sum += dis2[i];
      }
    mean = sum / N2;
    sum = 0;
    for( int i = 0; i < N2; i++ )
      {
      sum += vcl_pow(dis2[i] - mean, 2);
      }
    sd = vcl_sqrt( sum / ( N2 - 1 ) );
    upperT = mean + 3 * sd;
    std::cout << "Mis-match threshold2: " << upperT << " Mean " << mean << std::endl;
    index = 0;
    while( !fiberList2.empty() )
      {
      if( dis2[index] < lowerT )
        {
        merC++;
        }
      else if( dis2[index] < upperT )
        {
        mergedFiberList.push_back( fiberList2.front() );
        }
      else
        {
        disC++;
        }
      fiberList2.pop_front();
      index++;
      }

    std::cout << "Done!" << std::endl;
    finish = clock();
    float duration = (float)( finish - start ) / CLOCKS_PER_SEC;
    std::cout << "Time: " << duration << " seconds" << std::endl;
    std::cout << std::endl;
    std::cout << "original fiber: " << N1 << std::endl;
    std::cout << "fiber added:    " << N2 << std::endl;
    std::cout << "fibers merged:  " << merC << std::endl;
    std::cout << "fibers discard: " << disC << std::endl;
    std::cout << "Total fibers:   " << mergedFiberList.size() << std::endl;
    }

  return mergedFiberList;
}

void MinimumDistanceBetweenFiberGroups(VTKFiberListType fiberList1,
                                       VTKFiberListType fiberList2,
                                       TVector & dis)
{
  typedef itk::IdentityTransform<double> TransformType;
  TransformType::Pointer transform = TransformType::New();
  typedef itk::EuclideanDistancePointMetric<PointSetType, PointSetType> DistanceType;
  DistanceType::TransformParametersType parameter( transform->GetNumberOfParameters() );
  parameter.fill(0);

  VTKFiberListType::iterator fiberIt1, fiberIt2;
  int                        index = 0;
  for( fiberIt1 = fiberList1.begin(); fiberIt1 != fiberList1.end(); index++, fiberIt1++ )
    {
    vtkPolyData *      fiber1 = *fiberIt1;
    std::vector<float> avrDis;
    avrDis.clear();
    for( fiberIt2 = fiberList2.begin(); fiberIt2 != fiberList2.end(); fiberIt2++ )
      {
      vtkPolyData *         fiber2 = *fiberIt2;
      PointSetType::Pointer pSet1 = PolyDataToPointSet(fiber1);
      PointSetType::Pointer pSet2 = PolyDataToPointSet(fiber2);
      DistanceType::Pointer distance = DistanceType::New();
      distance->SetMovingPointSet(pSet1);
      distance->SetFixedPointSet(pSet2);
      distance->SetTransform(transform);
      DistanceType::MeasureType value = distance->GetValue(parameter);
      int                       n = distance->GetNumberOfValues();
      float                     avr = 0;
      for( int i = 0; i < n; i++ )
        {
        avr += value[i];
        }
      avr /= n;
      avrDis.push_back(avr);
      }
    float min = 999;
    for( unsigned int i = 0; i < avrDis.size(); i++ )
      {
      if( min > avrDis[i] )
        {
        min = avrDis[i];
        }
      }
    dis[index] = min;
    }
}

void SaveFiber(const char *fileName, VTKFiberListType & fiberList)
{
  std::cout << std::endl;
  std::cout << "Saving fibers...." << std::endl;
  std::cout << fileName << std::endl;
  vtkAppendPolyData *        append = vtkAppendPolyData::New();
  VTKFiberListType::iterator fiberIt;

  for( fiberIt = fiberList.begin(); fiberIt != fiberList.end(); fiberIt++ )
    {
    append->AddInput(*fiberIt);
    }

  vtkXMLPolyDataWriter *writer = vtkXMLPolyDataWriter::New();
  writer->SetFileName(fileName);
  writer->SetInput( append->GetOutput() );
  writer->SetDataModeToAppended();
  writer->EncodeAppendedDataOff();
  writer->SetDataModeToAscii();
  writer->SetCompressor(0);
  try
    {
    writer->Write();
    std::cout << "Done!" << std::endl;
    std::cout << "Number of fibers saved: " << append->GetOutput()->GetNumberOfLines() << std::endl;
    }
  catch( ... )
    {
    std::cout << "Error while saving fiber file." << std::endl;
    throw;
    }
}

PointSetType::Pointer PolyDataToPointSet(vtkPolyData *fiber)
{
  const int               npts = fiber->GetNumberOfPoints();
  PointSetType::Pointer   pSet = PointSetType::New();
  vtkFloatingPointType    p[3];
  PointSetType::PointType point;

  for( int i = 0; i < npts; i++ )
    {
    fiber->GetPoint(i, p);
    point[0] = p[0]; point[1] = p[1]; point[2] = p[2];
    pSet->SetPoint(i, point);
    }
  return pSet;
}

void OpenFiber(const char *fileName, VTKFiberListType & fiberList, bool reverse)
{
  if( fileName )
    {
    vtkXMLPolyDataReader *reader = vtkXMLPolyDataReader::New();
    reader->SetFileName(fileName);
    try
      {
      std::cout << std::endl;
      std::cout << "Loading fibers...." << std::endl;
      std::cout << fileName << std::endl;
      reader->Update();
      fiberList.clear();
      vtkPolyData * data  = reader->GetOutput();
      vtkCellArray *lines = data->GetLines();
      vtkDataArray *ai    = data->GetPointData()->GetScalars("Anisotropy");
      vtkDataArray *curv  = data->GetPointData()->GetScalars("Curvature");
      vtkIdType     npts, *pts;
      for( lines->InitTraversal(); lines->GetNextCell(npts, pts); )
        {
        vtkPoints *    points = vtkPoints::New();
        vtkFloatArray *scalars1 = vtkFloatArray::New();
        scalars1->SetName("Anisotropy");
        vtkFloatArray *scalars2 = vtkFloatArray::New();
        scalars2->SetName("Curvature");
        vtkCellArray *cells  = vtkCellArray::New();
        vtkPolyData * pdata  = vtkPolyData::New();
        cells->InsertNextCell(npts);
        for( int i = 0; i < npts; i++ )
          {
          if( reverse )
            {
            points->InsertPoint( i, data->GetPoint(pts[npts - 1 - i]) );
            if( ai )
              {
              scalars1->InsertTuple1( i, ai->GetTuple1(pts[npts - 1 - i]) );
              }
            if( curv )
              {
              scalars2->InsertTuple1( i, curv->GetTuple1(pts[npts - 1 - i]) );
              }
            }
          else
            {
            points->InsertPoint( i, data->GetPoint(pts[i]) );
            if( ai )
              {
              scalars1->InsertTuple1( i, ai->GetTuple1(pts[i]) );
              }
            if( curv )
              {
              scalars2->InsertTuple1( i, curv->GetTuple1(pts[i]) );
              }
            }
          cells->InsertCellPoint(i);
          }
        pdata->SetPoints(points);
        pdata->SetLines(cells);
        if( ai )
          {
          pdata->GetPointData()->AddArray(scalars1);
          pdata->GetPointData()->SetActiveScalars("Anisotropy");
          }
        if( curv )
          {
          pdata->GetPointData()->AddArray(scalars2);
          }
        fiberList.push_back(pdata);
        }
      std::cout << "Done!" << std::endl;
      std::cout << "Number of fibers loaded: " << fiberList.size() << std::endl;
      }
    catch( ... )
      {
      std::cout << "Error while reading fiber file." << std::endl;
      throw;
      }
    }
}
