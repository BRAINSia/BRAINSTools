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

int main(int argc, char *argv[])
{
  PARSE_ARGS;

  const int dimension = 3;
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
