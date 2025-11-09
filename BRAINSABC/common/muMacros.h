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
#ifndef __muMacros_h
#define __muMacros_h

#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#include <sstream>
#include <string>

#define muEcho(varname) std::cout << #varname << " = " << (varname) << std::flush << std::endl;

#define muStringMacro(strname, s) \
  std::string strname;            \
  {                               \
    std::ostringstream outss;     \
    outss << "" s << std::ends;   \
    (strname) = outss.str();      \
  }

#define muSelfFilterMacro(filter, obj) \
  {                                    \
    (filter)->SetInput(obj);           \
    iterator copy                      \
  }

#define muReadMacro(type, filename, image)                   \
  {                                                          \
    using ReaderType = itk::ImageFileReader<type>;           \
    typename ReaderType::Pointer reader = ReaderType::New(); \
    reader->SetFileName(filename);                           \
    reader->Update();                                        \
    (image) = reader->GetOutput();                           \
  }

#define muWriteMacro(type, filename, image)                  \
  {                                                          \
    using WriterType = itk::ImageFileWriter<type>;           \
    typename WriterType::Pointer writer = WriterType::New(); \
    writer->UseCompressionOn();                              \
    writer->SetFileName(filename);                           \
    writer->SetInput(image);                                 \
    writer->Update();                                        \
  }

#endif
