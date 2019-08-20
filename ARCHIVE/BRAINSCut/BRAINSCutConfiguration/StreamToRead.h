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
#include "itkMacro.h" //Needed for override

using FileStreamType = std::fstream;

class StreamToRead : public FileToRead< FileStreamType * >
{
public:
  ITK_DISALLOW_COPY_AND_ASSIGN( StreamToRead );

  using Self = StreamToRead;
  using Superclass = FileToRead< std::fstream * >;
  using OutputType = Superclass::OutputType;

  OutputType
  GetValue() const override;

  void
  Close() override;

  int
  PrintSelf( std::ostream & os, int indent ) const override;

  StreamToRead( const std::string & name, const std::string & filename );
  ~StreamToRead() override;

protected:
  StreamToRead() ITK_DELETED_FUNCTION;

private:
  OutputType m_F;
};
#endif // StreamToRead_h
