// Simple utility functions

#ifndef __muUtils_h
#define __muUtils_h

#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#include "muException.h"

#include <iostream>

#define muAlloc(n, T) muAlloc_func<T>(n, __FILE__, __LINE__)

template <class T>
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

template <class TImage>
typename TImage::Pointer
readImage(const char *fn)
{
  typedef itk::ImageFileReader<TImage> ReaderType;
  typename ReaderType::Pointer reader = ReaderType::New();

  reader->SetFileName(fn);
  reader->Update();

  return reader->GetOutput();
}

template <class TImage>
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

// template <class TImage>
// void
// writeImageAsByte
// {
//  typedef itk::Cast
//  typedef itk::ImageFileWriter<TImage> WriterType;
// }

// void
// writeImageAsShort

#endif
