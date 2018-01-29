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

 Copyright (c) University of Iowa Department of Radiology. All rights reserved.
 See GTRACT-Copyright.txt or http://mri.radiology.uiowa.edu/copyright/GTRACT-Copyright.txt
 for details.

 This software is distributed WITHOUT ANY WARRANTY; without even
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE.  See the above copyright notices for more information.

=========================================================================*/

#include "vtksys/SystemTools.hxx"
#include "vtkStructuredGridReader.h"
#include "vtkXMLStructuredGridReader.h"
#include "vtkTalairachConversion.h"
#include "BRAINSTalairachMaskCLP.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include <BRAINSCommonLib.h>

int main(int argc, char *argv[])
{
  PARSE_ARGS;
  BRAINSRegisterAlternateIO();

  constexpr int dimension = 3;
  typedef itk::Image<unsigned char, dimension> ImageType;
  typedef itk::ImageFileReader<ImageType>      ImageReaderType;

  ImageReaderType::Pointer reader = ImageReaderType::New();
  reader->SetFileName(inputVolume);
  reader->Update();

  vtkStructuredGrid *talairach;
  const std::string extension = vtksys::SystemTools::LowerCase( vtksys::SystemTools::GetFilenameLastExtension( talairachParameters ) );

  if( extension == ".vtk" )
    {
    vtkStructuredGridReader *gridReader = vtkStructuredGridReader::New();
    gridReader->SetFileName( talairachParameters.c_str() );
    gridReader->Update();
    talairach = gridReader->GetOutput();
    }
  else
    {
    vtkXMLStructuredGridReader *gridReader = vtkXMLStructuredGridReader::New();
    gridReader->SetFileName( talairachParameters.c_str() );
    gridReader->Update();
    talairach = gridReader->GetOutput();
    }

  vtkTalairachConversion *tConv = vtkTalairachConversion::New();
  tConv->SetImageInformation( reader->GetOutput() );
  tConv->SetTalairachGrid( talairach );
  if( hemisphereMode == "right" )
    {
    tConv->SetHemisphereMode( vtkTalairachConversion::right );
    }
  else if( hemisphereMode == "left" )
    {
    tConv->SetHemisphereMode( vtkTalairachConversion::left );
    }
  else
    {
    tConv->SetHemisphereMode( vtkTalairachConversion::both );
    }
  tConv->SetSegmentationMode( expand );

  ifstream    fin( talairachBox.c_str() );
  std::string line;
  getline(fin, line);

  while( !line.empty() )
    {
    tConv->AddTalairachBox( line );
    line.clear();
    getline(fin, line);
    }

  tConv->Update();

  typedef itk::ImageFileWriter<ImageType> ImageWriterType;
  ImageWriterType::Pointer writer = ImageWriterType::New();
  writer->SetFileName( outputVolume );
  writer->SetInput( tConv->GetImage() );
  writer->Update();

  return EXIT_SUCCESS;
}
