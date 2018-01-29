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
#include <iostream>
#include <fstream>
#include <sstream>

#include "itkIO.h"
#include <BRAINSCommonLib.h>

#include "DumpBinaryTrainingVectorsCLP.h"

#define MAX_LINE_SIZE 1000

static const float LineGuard = 1234567.0;

void
ReadHeader(const char *fname,
           unsigned int & InputVectorSize,
           unsigned int & OutputVectorSize,
           unsigned int & NumberTrainingVectorsFromFile)
{
  std::ifstream filestr;

  filestr.open(fname);
  if( !filestr.is_open() )
    {
    itkGenericExceptionMacro(<< "Error: Could not open ANN vector file::" << fname );
    }
  for( int tags = 0; tags < 3; tags++ )
    {
    std::string temp;
    char        currentline[MAX_LINE_SIZE];
    filestr.getline(currentline, MAX_LINE_SIZE - 1);
    std::istringstream iss(currentline, std::istringstream::in);
    iss >> temp;
    if( temp == "IVS" )
      {
      iss >> InputVectorSize;
      }
    else if( temp == "OVS" )
      {
      iss >> OutputVectorSize;
      }
    else if( temp == "TVC" )
      {
      iss >> NumberTrainingVectorsFromFile;
      }
    }
  filestr.close();
}

int
main(int argc, char * *argv)
{
  PARSE_ARGS;
  BRAINSRegisterAlternateIO();

  unsigned int InputVectorSize;
  unsigned int OutputVectorSize;
  unsigned int NumberTrainingVectorsFromFile;
  int          result(0);

  if( inputVectorFilename.empty() ||
      inputHeaderFilename.empty() )
    {
    std::cout << "Error: Required inputs are missing:  --inputVectorFilename and/or --inputHeaderFilename "
              << std::endl;
    std::exit(EXIT_FAILURE);
    }
  // read header
  try
    {
    ReadHeader( inputHeaderFilename.c_str(),
                InputVectorSize,
                OutputVectorSize,
                NumberTrainingVectorsFromFile);
    std::cout << "Input Vector  Size "  << InputVectorSize
              << "Output Vector Size "  << OutputVectorSize
              << " #Training Vectors " << NumberTrainingVectorsFromFile
              << std::endl;
    }
  catch( std::exception & ex )
    {
    std::cerr << ex.what() << std::endl;
    return EXIT_FAILURE;
    }
  // read binary
  // constexpr unsigned int OutputVectorSize = 2;
  constexpr unsigned int SentinalValueSize = 1;
  std::ifstream      binfile;
  binfile.open(inputVectorFilename.c_str(), std::ios::in | std::ios::binary);
  if( !binfile.is_open() )
    {
    std::cerr << "Can't open " << inputVectorFilename << std::endl;
    return EXIT_FAILURE;
    }
  unsigned int recordsize =
    ( InputVectorSize + OutputVectorSize
      + SentinalValueSize ) * sizeof( float );
  float *buf =
    new float[InputVectorSize + OutputVectorSize + SentinalValueSize];
  for( unsigned int j = 0; j < OutputVectorSize; j++ )
    {
    std::cout << " Out" << ( j + 1 );
    }
  for( unsigned int j = 0; j < InputVectorSize; j++ )
    {
    std::cout << " In" << ( j + 1 );
    }
  std::cout << std::endl;
  for( unsigned int i = 0; i < NumberTrainingVectorsFromFile; i++ )
    {
    if( binfile.eof() )
      {
      std::cerr << "Premature end of file at record "
                << i << std::endl;
      result = 1;
      break;
      }
    binfile.read( (char *)buf, recordsize );
    for( unsigned int j = 0;
         j <= OutputVectorSize + InputVectorSize;
         j++ )
      {
      std::cout << buf[j] << " ";
      }
    std::cout << std::endl;
    if( buf[OutputVectorSize + InputVectorSize + SentinalValueSize - 1] !=
        LineGuard )
      {
      std::cerr << "Record not properly terminated by sentinel value :::: "
                << buf[OutputVectorSize + InputVectorSize + SentinalValueSize - 1]
                << " != " << LineGuard
                << std::endl;
      result = 1;
      break;
      }
    }
  delete[] buf;
  return result;
}
