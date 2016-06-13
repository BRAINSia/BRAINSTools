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
#ifndef StreamToRead_h
#define StreamToRead_h
#include "FileToRead.h"
#include <iostream>
#include <fstream>
#include "itkMacro.h" //Needed for ITK_OVERRIDE

typedef std::fstream FileStreamType;

class StreamToRead :
  public FileToRead<FileStreamType *>
{
public:
  typedef StreamToRead               Self;
  typedef FileToRead<std::fstream *> Superclass;
  typedef Superclass::OutputType     OutputType;

  virtual OutputType GetValue() const ITK_OVERRIDE;

  virtual void Close() ITK_OVERRIDE;

  virtual int PrintSelf(std::ostream & os, int indent) const ITK_OVERRIDE;

  StreamToRead(const std::string & name, const std::string & filename);
  virtual ~StreamToRead();
protected:
  StreamToRead() ITK_DELETED_FUNCTION;
  ITK_DISALLOW_COPY_AND_ASSIGN(StreamToRead);

private:
  OutputType m_F;
};
#endif // StreamToRead_h
