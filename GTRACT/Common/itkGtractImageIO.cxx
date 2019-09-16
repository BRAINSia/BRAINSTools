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

#ifndef __itkGtractImageIO_cxx
#define __itkGtractImageIO_cxx

#include "itkGtractImageIO.h"

#include <iostream>

namespace itk
{
GtractImageIO ::GtractImageIO() {}

void
GtractImageIO::SetDicomDirectory(char * dicomDir)
{
  m_DicomDirectory = dicomDir;
}

void
GtractImageIO::SetDicomDirectory(std::string dicomDir)
{
  m_DicomDirectory = dicomDir;
}

void
GtractImageIO::SetFileName(char * fileName)
{
  m_FileName = fileName;
}

void
GtractImageIO::SetFileName(std::string fileName)
{
  m_FileName = fileName;
}

void
GtractImageIO::SetDicomSeriesUID(char * UID)
{
  m_DicomSeriesUID = UID;
}

void
GtractImageIO::SetDicomSeriesUID(std::string UID)
{
  m_DicomSeriesUID = UID;
}

void
GtractImageIO::Load3dShortImage()
{
  using FileReaderType = itk::ImageFileReader<Short3dImageType>;
  FileReaderType::Pointer reader = FileReaderType::New();
  std::cout << "Loading image " << m_FileName << " ...." << std::endl;
  reader->SetFileName(m_FileName.c_str());
  reader->Update();

  m_Short3dImage = reader->GetOutput();
}

void
GtractImageIO::Load4dShortImage()
{
  using FileReaderType = itk::ImageFileReader<Short4dImageType>;
  FileReaderType::Pointer reader = FileReaderType::New();
  std::cout << "Loading image " << m_FileName << " ...." << std::endl;
  reader->SetFileName(m_FileName.c_str());
  reader->Update();

  m_Short4dImage = reader->GetOutput();
}

void
GtractImageIO::Load3dFloatImage()
{
  using FileReaderType = itk::ImageFileReader<Float3dImageType>;
  FileReaderType::Pointer reader = FileReaderType::New();
  std::cout << "Loading image " << m_FileName << " ...." << std::endl;
  reader->SetFileName(m_FileName.c_str());
  reader->Update();

  m_Float3dImage = reader->GetOutput();
}

void
GtractImageIO::Load3dRgbImage()
{
  using FileReaderType = itk::ImageFileReader<Rgb3dImageType>;
  FileReaderType::Pointer reader = FileReaderType::New();
  std::cout << "Loading image " << m_FileName << " ...." << std::endl;
  reader->SetFileName(m_FileName.c_str());
  reader->Update();

  m_Rgb3dImage = reader->GetOutput();
}

void
GtractImageIO::LoadTensorImage()
{
  using FileReaderType = itk::ImageFileReader<TensorImageType>;
  FileReaderType::Pointer reader = FileReaderType::New();
  std::cout << "Loading image " << m_FileName << " ...." << std::endl;
  reader->SetFileName(m_FileName.c_str());
  reader->Update();

  m_TensorImage = reader->GetOutput();
}

void
GtractImageIO::Save3dShortImage()
{
  using FileWriterType = itk::ImageFileWriter<Short3dImageType>;
  FileWriterType::Pointer writer = FileWriterType::New();
  writer->UseCompressionOn();
  std::cout << "Saving image " << m_FileName << " ...." << std::endl;
  writer->SetInput(m_Short3dImage);
  writer->SetFileName(m_FileName.c_str());
  writer->Update();
  std::cout << "Done!" << std::endl;
}

void
GtractImageIO::Save4dShortImage()
{
  using FileWriterType = itk::ImageFileWriter<Short4dImageType>;
  FileWriterType::Pointer writer = FileWriterType::New();
  writer->UseCompressionOn();
  std::cout << "Saving image " << m_FileName << " ...." << std::endl;
  writer->SetInput(m_Short4dImage);
  writer->SetFileName(m_FileName.c_str());
  writer->Update();
  std::cout << "Done!" << std::endl;
}

void
GtractImageIO::Save3dFloatImage()
{
  using FileWriterType = itk::ImageFileWriter<Float3dImageType>;
  FileWriterType::Pointer writer = FileWriterType::New();
  writer->UseCompressionOn();
  std::cout << "Saving image " << m_FileName << " ...." << std::endl;
  writer->SetInput(m_Float3dImage);
  writer->SetFileName(m_FileName.c_str());
  writer->Update();
  std::cout << "Done!" << std::endl;
}

void
GtractImageIO::Save3dRgbImage()
{
  using FileWriterType = itk::ImageFileWriter<Rgb3dImageType>;
  FileWriterType::Pointer writer = FileWriterType::New();
  writer->UseCompressionOn();
  std::cout << "Saving image " << m_FileName << " ...." << std::endl;
  writer->SetInput(m_Rgb3dImage);
  writer->SetFileName(m_FileName.c_str());
  writer->Update();
  std::cout << "Done!" << std::endl;
}

void
GtractImageIO::SaveTensorImage()
{
  using FileWriterType = itk::ImageFileWriter<TensorImageType>;
  FileWriterType::Pointer writer = FileWriterType::New();
  writer->UseCompressionOn();
  std::cout << "Saving image " << m_FileName << " ...." << std::endl;
  writer->SetInput(m_TensorImage);
  writer->SetFileName(m_FileName.c_str());
  writer->Update();
  std::cout << "Done!" << std::endl;
}

void
GtractImageIO::Load3dDICOMSeries()
{
  using ReaderType = itk::ImageSeriesReader<Short3dImageType>;
  ReaderType::Pointer       reader = ReaderType::New();
  itk::GDCMImageIO::Pointer dicomIO = itk::GDCMImageIO::New();

  /* Generate a list of Series UIDs */
  itk::GDCMSeriesFileNames::Pointer FileNameGenerator;
  FileNameGenerator = itk::GDCMSeriesFileNames::New();
  FileNameGenerator->SetUseSeriesDetails(false);
  FileNameGenerator->SetDirectory(m_DicomDirectory.c_str());

  reader->SetFileNames(FileNameGenerator->GetFileNames(m_DicomSeriesUID.c_str()));
  reader->SetImageIO(dicomIO);

  std::cout << "Loading dicom images...." << std::endl;
  reader->Update();
  std::cout << "Done!" << std::endl;

  m_Short3dImage = reader->GetOutput();
}

void
GtractImageIO::Load4dDICOMSeries()
{
  using ReaderType = itk::ImageSeriesReader<Short4dImageType>;
  ReaderType::Pointer       reader = ReaderType::New();
  itk::GDCMImageIO::Pointer dicomIO = itk::GDCMImageIO::New();

  /* Generate a list of Series UIDs */
  itk::GDCMSeriesFileNames::Pointer FileNameGenerator;
  FileNameGenerator = itk::GDCMSeriesFileNames::New();
  FileNameGenerator->SetUseSeriesDetails(true);
  FileNameGenerator->SetDirectory(m_DicomDirectory.c_str());

  reader->SetFileNames(FileNameGenerator->GetFileNames(m_DicomSeriesUID.c_str()));
  reader->SetImageIO(dicomIO);

  std::cout << "Loading dicom images...." << std::endl;
  reader->Update();
  std::cout << "Done!" << std::endl;

  m_Short4dImage = reader->GetOutput();
}
} // end namespace itk
#endif
