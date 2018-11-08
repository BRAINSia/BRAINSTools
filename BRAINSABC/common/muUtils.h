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
// Simple utility functions

#ifndef __muUtils_h
#define __muUtils_h

#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#include "muException.h"

#include <iostream>

#define muAlloc(n, T) muAlloc_func<T>(n, __FILE__, __LINE__)

template <typename T>
T *
muAlloc_func(unsigned int n, const char *s, int line)
{
  T *array = new T[n];

  if( array == NULL )
    {
    muExceptionMacro(<< "muAlloc: Failed at " << s << ": " << line);
    }
  return array;
}

template <typename TImage>
typename TImage::Pointer
readImage(const char *fn)
{
  typedef itk::ImageFileReader<TImage> ReaderType;
  typename ReaderType::Pointer reader = ReaderType::New();

  reader->SetFileName(fn);
  reader->Update();

  return reader->GetOutput();
}

template <typename TImage>
void
writeImage(const char *fn, const TImage *ip)
{
  typedef itk::ImageFileWriter<TImage> WriterType;
  typename WriterType::Pointer writer = WriterType::New();
  writer->UseCompressionOn();

  writer->SetFileName(fn);
  writer->SetInput(ip);

  writer->Update();
}

// template <typename TImage>
// void
// writeImageAsByte
// {
//  typedef itk::Cast
//  typedef itk::ImageFileWriter<TImage> WriterType;
// }

// void
// writeImageAsShort

#endif
