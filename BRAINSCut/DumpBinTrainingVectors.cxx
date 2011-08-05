#include "Utilities.h"
#include <iostream>
#include <fstream>
#include <sstream>

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
    std::cout << "Error: Could not open ANN vector file" << std::endl;
    exit(1);
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
  unsigned int InputVectorSize;
  unsigned int OutputVectorSize;
  unsigned int NumberTrainingVectorsFromFile;
  int          result(0);

  if( argc < 3 )
    {
    std::cerr << "DumpBinTrainingVectors hdr-name bin-name [record-offset]"
              << std::endl;
    exit(1);
    }
  // read header
  ReadHeader(argv[1],
             InputVectorSize,
             OutputVectorSize,
             NumberTrainingVectorsFromFile);
  std::cout << "Input Vector  Size "  << InputVectorSize
            << "Output Vector Size "  << OutputVectorSize
            << " #Training Vectors " << NumberTrainingVectorsFromFile
            << std::endl;
  // read binary
  // const unsigned int OutputVectorSize=2;
  const unsigned int SentinalValueSize = 1;
  std::ifstream      binfile;
  binfile.open(argv[2], std::ios::in | std::ios::binary);
  if( !binfile.is_open() )
    {
    std::cerr << "Can't open " << argv[2] << std::endl;
    exit(1);
    }
  unsigned skip = 0;
  if( argc > 3 )
    {
    skip = atoi(argv[3]);
    }
  if( skip > 0 )
    {
    std::cout << "Skipping " << skip << " vectors"
              << std::endl;
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
    for( unsigned int j = 0; j < OutputVectorSize; j++ )
      {
      std::cout << buf[j] << " ";
      }
    binfile.read( (char *)buf, recordsize );
    for( unsigned int j = OutputVectorSize;
         j < OutputVectorSize + InputVectorSize;
         j++ )
      {
      std::cout << buf[j] << " ";
      }

    std::cout << std::endl;
    if( buf[OutputVectorSize + InputVectorSize + SentinalValueSize - 1] !=
        AUTOSEG_VEC_SENTINEL )
      {
      std::cerr << "Record not properly terminated by sentinel value"
                << std::endl;
      result = 1;
      break;
      }
    }
  delete[] buf;
  exit(result);
}
