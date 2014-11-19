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

#include "BRAINSTalairachCLP.h"
#include "vtkTalairachGrid.h"
#include <string>
#include "vtkXMLStructuredGridWriter.h"
#include "vtkStructuredGridWriter.h"
#include "vtksys/SystemTools.hxx"
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkPoint.h"
#include "itkImageIteratorWithIndex.h"
#include <BRAINSCommonLib.h>

int main(int argc, char *argv[])
{
  PARSE_ARGS;
  BRAINSRegisterAlternateIO();

  if( ( AC.size() != 3 ) || ( PC.size() != 3 ) || ( IRP.size() != 3 ) || ( SLA.size() != 3 ) )
    {
    std::cout << "Error: The AC, PC, IRP, and SLA points must all have three values" << std::endl;
    return 1;
    }

  const bool debug = true;
  if( debug )
    {
    std::cout << "=====================================================" << std::endl;
    std::cout << "AC: " <<   AC[0] << " " <<  AC[1] << " " <<  AC[2] << std::endl;
    std::cout << "PC: " <<   PC[0] << " " <<  PC[1] << " " <<  PC[2] << std::endl;
    std::cout << "IRP: " << IRP[0] << " " << IRP[1] << " " << IRP[2] << std::endl;
    std::cout << "SLA: " <<  SLA[0] << " " << SLA[1] << " " << SLA[2] << std::endl;
    std::cout << "ACisIndex: " <<  ACisIndex << std::endl;
    std::cout << "PCisIndex: " <<  PCisIndex << std::endl;
    std::cout << "IRPisIndex: " <<  IRPisIndex << std::endl;
    std::cout << "SLAisIndex: " <<  SLAisIndex << std::endl;
    std::cout << "Image: " <<  inputVolume << std::endl;
    std::cout << "Talairach Grid: " <<  outputGrid << std::endl;
    std::cout << "Talairach Box: " <<  outputBox << std::endl;
    std::cout << "=====================================================" << std::endl;
    }

  const int dimension = 3;
  typedef itk::ContinuousIndex<double, 3>      ContinuousIndexType;
  typedef itk::Image<unsigned char, dimension> ImageType;
  typedef itk::ImageFileReader<ImageType>      ImageReaderType;
  ImageReaderType::Pointer reader = ImageReaderType::New();
  reader->SetFileName(inputVolume);
  reader->Update();

  double ACpoint[3];
  double PCpoint[3];
  double IRPpoint[3];
  double SLApoint[3];
  for( int i = 0; i < 3; i++ )
    {
    ACpoint[i] = AC[i];
    PCpoint[i] = PC[i];
    IRPpoint[i] = IRP[i];
    SLApoint[i] = SLA[i];
    }

  vtkTalairachGrid *tGrid = vtkTalairachGrid::New();

  if( ACisIndex )
    {
    ImageType::PointType point;
    ContinuousIndexType  pixelIndex;
    pixelIndex[0] = ACpoint[0];
    pixelIndex[1] = ACpoint[1];
    pixelIndex[2] = ACpoint[2];
    ( reader->GetOutput() )->TransformContinuousIndexToPhysicalPoint(pixelIndex, point);
    ACpoint[0] = point[0];
    ACpoint[1] = point[1];
    ACpoint[2] = point[2];
    }

  if( PCisIndex )
    {
    ImageType::PointType point;
    ContinuousIndexType  pixelIndex;
    pixelIndex[0] = PCpoint[0];
    pixelIndex[1] = PCpoint[1];
    pixelIndex[2] = PCpoint[2];
    ( reader->GetOutput() )->TransformContinuousIndexToPhysicalPoint(pixelIndex, point);
    PCpoint[0] = point[0];
    PCpoint[1] = point[1];
    PCpoint[2] = point[2];
    }

  if( IRPisIndex )
    {
    ImageType::PointType point;
    ContinuousIndexType  pixelIndex;
    pixelIndex[0] = IRPpoint[0];
    pixelIndex[1] = IRPpoint[1];
    pixelIndex[2] = IRPpoint[2];
    ( reader->GetOutput() )->TransformContinuousIndexToPhysicalPoint(pixelIndex, point);
    std::cout << "Mapped IRP: " << point << std::endl;
    IRPpoint[0] = point[0];
    IRPpoint[1] = point[1];
    IRPpoint[2] = point[2];
    }

  if( SLAisIndex )
    {
    ImageType::PointType point;
    ContinuousIndexType  pixelIndex;
    pixelIndex[0] = SLApoint[0];
    pixelIndex[1] = SLApoint[1];
    pixelIndex[2] = SLApoint[2];
    ( reader->GetOutput() )->TransformContinuousIndexToPhysicalPoint(pixelIndex, point);
    std::cout << "Mapped SLA: " << point << std::endl;
    SLApoint[0] = point[0];
    SLApoint[1] = point[1];
    SLApoint[2] = point[2];
    }

  tGrid->SetACPoint(ACpoint);
  tGrid->SetPCPoint(PCpoint);

  if( IRPpoint[0] > SLApoint[0] )
    {
    double tmp;
    tmp = SLApoint[0];
    SLApoint[0] = IRPpoint[0];
    IRPpoint[0] = tmp;
    }

  if( IRPpoint[1] > SLApoint[1] )
    {
    double tmp;
    tmp = SLApoint[1];
    SLApoint[1] = IRPpoint[1];
    IRPpoint[1] = tmp;
    }

  if( IRPpoint[2] > SLApoint[2] )
    {
    double tmp;
    tmp = SLApoint[2];
    SLApoint[2] = IRPpoint[2];
    IRPpoint[2] = tmp;
    }

  tGrid->SetIRPPoint(IRPpoint);
  tGrid->SetSLAPoint(SLApoint);
  tGrid->EstablishBoundingBoxGrid();
  tGrid->EstablishTalairachGrid();
  tGrid->Update();

  vtkStructuredGrid *talairach;

  if( outputGrid.length() > 0 )
    {
    talairach = tGrid->GetTalairachGrid();

    std::string extension = vtksys::SystemTools::LowerCase(
        vtksys::SystemTools::GetFilenameLastExtension( outputGrid ) );

    if( extension == ".vtk" )
      {
      vtkStructuredGridWriter *writer = vtkStructuredGridWriter::New();
#if (VTK_MAJOR_VERSION < 6)
      writer->SetInput( talairach );
#else
      writer->SetInputData( talairach );
#endif
      writer->SetFileName( outputGrid.c_str() );
      writer->Update();
      }
    else
      {
      vtkXMLStructuredGridWriter *writer = vtkXMLStructuredGridWriter::New();
#if (VTK_MAJOR_VERSION < 6)
      writer->SetInput( talairach );
#else
      writer->SetInputData( talairach );
#endif
      writer->SetFileName( outputGrid.c_str() );
      writer->Update();
      }
    }

  if( outputBox.length() > 0 )
    {
    talairach = tGrid->GetBoundingBoxGrid();

    std::string extension = vtksys::SystemTools::LowerCase(
        vtksys::SystemTools::GetFilenameLastExtension( outputBox ) );

    if( extension == ".vtk" )
      {
      vtkStructuredGridWriter *writer = vtkStructuredGridWriter::New();
#if (VTK_MAJOR_VERSION < 6)
      writer->SetInput( talairach );
#else
      writer->SetInputData( talairach );
#endif
      writer->SetFileName( outputBox.c_str() );
      writer->Update();
      }
    else
      {
      vtkXMLStructuredGridWriter *writer = vtkXMLStructuredGridWriter::New();
#if (VTK_MAJOR_VERSION < 6)
      writer->SetInput( talairach );
#else
      writer->SetInputData( talairach );
#endif
      writer->SetFileName( outputBox.c_str() );
      writer->Update();
      }
    }

  return 0;
}
