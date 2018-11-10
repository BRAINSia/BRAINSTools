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
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImage.h"
#include "itkPoint.h"
#include "itkImageRegionConstIterator.h"

#include "vtkVersion.h"
#include "vtkPolyDataReader.h"
#include "vtkPolyDataWriter.h"
#include "vtkPolyData.h"
#include "vtkPointData.h"
#include "vtkPoints.h"
#include "vtkFloatArray.h"
#include "vtkTransform.h"
#include "vtkTransformPolyDataFilter.h"
#include "vtkSmartPointer.h"

#include "SurfaceColorCLP.h"
#include <BRAINSCommonLib.h>

int main( int argc, char * argv[] )
{
  PARSE_ARGS;
  BRAINSRegisterAlternateIO();

  std::cout << "---------------------------------------------------" << std::endl;
  std::cout << "Apply " << numOfLabels << " labels from: " << std::endl;
  std::cout << labelMapFile << std::endl;
  std::cout << "to the surface of: " << std::endl;
  std::cout << inputSurfaceFile << std::endl;
  std::cout << "Search labels for surface points in radius of " << radius << std::endl;
  if( increaseRadius )
    {
    std::cout << "Increase the radius until label is found for every surface point" << std::endl;
    }
  else
    {
    std::cout << "Searching radius is fixed" << std::endl;
    }
  std::cout << "---------------------------------------------------" << std::endl;

  // read the surface
  vtkSmartPointer<vtkPolyDataReader> surfacereader = vtkSmartPointer<vtkPolyDataReader>::New();
  surfacereader->SetFileName(inputSurfaceFile.c_str() );
  surfacereader->Update();

  // transform the surface into ITK image coordinates
  vtkSmartPointer<vtkTransform> rasOrientation = vtkSmartPointer<vtkTransform>::New();
  const double                  orientationMatrix[16] = { -1,  0, 0, 0,
                                                          0, -1, 0, 0,
                                                          0,  0, 1, 0,
                                                          0,  0, 0, 0  };
  rasOrientation->SetMatrix( orientationMatrix );

  vtkSmartPointer<vtkTransformPolyDataFilter> resampleSurface = vtkSmartPointer<vtkTransformPolyDataFilter>::New();
#if (VTK_MAJOR_VERSION < 6)
  resampleSurface->SetInput( surfacereader->GetOutput() );
#else
  resampleSurface->SetInputData( surfacereader->GetOutput() );
#endif
  resampleSurface->SetTransform( rasOrientation );
  resampleSurface->Update();

  vtkPolyData *surface = resampleSurface->GetOutput();

  // define image type
  constexpr unsigned char dimension = 3;
  using PixelType = unsigned char;
  using ImageType = itk::Image<PixelType, dimension>;

  // read Inputimage
  using ImageReaderType = itk::ImageFileReader<ImageType>;
  ImageReaderType::Pointer imagereader = ImageReaderType::New();

  imagereader->SetFileName(labelMapFile.c_str() );
  imagereader->Update();

  ImageType::ConstPointer image = imagereader->GetOutput();

  ImageType::SizeType  size = image->GetBufferedRegion().GetSize();
  ImageType::IndexType start = image->GetBufferedRegion().GetIndex();

  vtkPoints *surfacepoint = surface->GetPoints();
  int        npoints = surfacepoint->GetNumberOfPoints();

  vtkSmartPointer<vtkFloatArray> label = vtkSmartPointer<vtkFloatArray>::New();
  double                         point[3];

  ImageType::IndexType pixelindex;
  using ConstIteratorType = itk::ImageRegionConstIterator<ImageType>;

  ImageType::RegionType            subregion;
  ImageType::RegionType::IndexType substart;
  ImageType::RegionType::SizeType  subsize;

  using PointType = itk::Point<double, ImageType::ImageDimension>;
  PointType impoint;
  for( int i = 0; i < npoints; i++ )
    {
    surface->GetPoint(i, point); // physical points

    impoint[0] = point[0];
    impoint[1] = point[1];
    impoint[2] = point[2];

    bool Isin = image->TransformPhysicalPointToIndex( impoint, pixelindex );

    if( !Isin )
      {
      std::cout << "There are points of input surface located outside of image." << std::endl;
      std::cout << "Quit." << std::endl;
      return EXIT_FAILURE;
      }

    ImageType::PixelType pixelValue = image->GetPixel( pixelindex );

    // instead of getting pixel directly
    // find the sizexsize neighborhood and take the most occurring non-zero value
    if( pixelValue == 0 )
      {
      bool Found = false;

      while( !Found )
        {
        substart[0] = pixelindex[0] - radius;
        substart[1] = pixelindex[1] - radius;
        substart[2] = pixelindex[2] - radius;

        if( substart[0] < start[0] )
          {
          substart[0] = start[0];
          }
        if( substart[1] < start[1] )
          {
          substart[1] = start[1];
          }
        if( substart[2] < start[2] )
          {
          substart[2] = start[2];
          }

        subsize[0] = radius * 2 + 1;
        subsize[1] = radius * 2 + 1;
        subsize[2] = radius * 2 + 1;

        subregion.SetIndex( substart );
        subregion.SetSize( subsize );

        ConstIteratorType it( image, subregion );

        int *counter = new int[numOfLabels];
        for( int ii = 0; ii < numOfLabels; ii++ )
          {
          counter[ii] = 0;
          }

        // save number of pixels with value i
        int newValue = 0, temp = 0;
        for( it.GoToBegin(); !it.IsAtEnd(); ++it )
          {
          ImageType::IndexType subIndex = it.GetIndex(); // no bounds checking

          if( (subIndex[0] >= start[0] && subIndex[0] <= start[0] + int(size[0]) - 1)
              && (subIndex[1] >= start[1] && subIndex[1] <= start[1] + int(size[1]) - 1)
              && (subIndex[2] >= start[2] && subIndex[2] <= start[2] + int(size[2]) - 1) )
            {
            ImageType::PixelType subValue = it.Get();
            if( subValue != 0 )
              {
              counter[subValue - 1] += 1;
              }
            }
          }
        for( int j = 0; j < numOfLabels; j++ )
          {
          if( counter[j] > temp )
            {
            temp = counter[j];
            newValue = j + 1;
            }
          }

        if( (newValue != 0) || (radius > 15.0) || (!increaseRadius) )
          {
          Found = true;
          }
        else if( newValue == 0 )
          {
          std::cout << "some point can not find valid label." << std::endl;
          std::cout << "Continue" << std::endl;
          radius += 1;
          std::cout << "radius: " << radius << std::endl;
          }

        pixelValue = newValue;
        delete [] counter;
        }
      }

    label->InsertValue( i, pixelValue );
    }

  label->SetName("LabelValue");
  if( surface->GetPointData()->GetScalars() == nullptr )
    {
    surface->GetPointData()->SetScalars(label);
    }
  else
    {
    surface->GetPointData()->AddArray(label);
    surface->GetPointData()->SetActiveScalars("LabelValue");
    }

  // Put the surface back into its original orientation
  vtkSmartPointer<vtkTransformPolyDataFilter> revertedSurface = vtkSmartPointer<vtkTransformPolyDataFilter>::New();
#if (VTK_MAJOR_VERSION < 6)
  revertedSurface->SetInput( surface );
#else
  revertedSurface->SetInputData( surface );
#endif
  revertedSurface->SetTransform( rasOrientation->GetInverse() );
  revertedSurface->Update();

  vtkSmartPointer<vtkPolyDataWriter> writer = vtkSmartPointer<vtkPolyDataWriter>::New();
#if (VTK_MAJOR_VERSION < 6)
  writer->SetInput(revertedSurface->GetOutput() );
#else
  writer->SetInputData(revertedSurface->GetOutput() );
#endif
  writer->SetFileName(outputSurfaceFile.c_str() );
  writer->SetFileTypeToASCII();
  writer->Update();

  return EXIT_SUCCESS;
}
