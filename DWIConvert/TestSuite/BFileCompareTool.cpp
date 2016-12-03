#include "stdio.h"
#include <iostream>
#include <vector>
#include <cmath>

typedef std::vector< std::vector< double > > DataVector;

static const int USAGE_FAILURE_CODE=1;
static const int FILE_READ_FAILURE_CODE=2;
static const int VALUE_TOLERANCE_FAILURE_CODE=4;

int readFileIntoBuffer(const char* filename, char** buffer){
  FILE* pFile = fopen(filename, "rb");
  if (NULL == pFile) {
    std::cerr<<"Sorry, can not open file: "<<filename<<std::endl;
    return FILE_READ_FAILURE_CODE;
  }

  fseek(pFile,0,SEEK_END);
  long nFileSize = ftell(pFile);
  fseek(pFile,0,SEEK_SET);

  *buffer = new char[nFileSize+1];
  long readSize = fread(*buffer, sizeof(char), (std::size_t) nFileSize, pFile);
  fclose(pFile);
  *(*buffer + readSize) = '\0';
  return EXIT_SUCCESS; //success
}

char* nextNoneSpacePointer(char* p){
  while(' ' == *p){
    ++p;
  }
  return p;
}

//notes: after using buffer, this function deletes buffer.
DataVector convertBufferIntoVector(char* buffer){
  DataVector ColumnData;
  char*p = nextNoneSpacePointer(buffer);
  while ('\0' != *p){
    std::vector< double > lineData;
    lineData.reserve(200);
    while ('\n' != *p && '\0' != *p){
      std::size_t dataOffset = 0;
      char *dataBeginning = p;
      while (' ' != *p && '\n' != *p && '\0' != *p) {
        ++dataOffset;
        ++p;
      }
      char *dataBuf = new char[dataOffset + 1];
      memcpy(dataBuf, dataBeginning, dataOffset);
      dataBuf[dataOffset] = '\0';
      double dataValue = atof(dataBuf);
      delete[] dataBuf;
      lineData.push_back(dataValue);
      p = nextNoneSpacePointer(p);
    }
    if (lineData.size() > 0) ColumnData.push_back(lineData);
    if ('\0' != *p) ++p;
  }
  delete[] buffer;
  return ColumnData;
}

//return: 0: all comparing data are close to each other
int compareRowColumn(const DataVector& data1, const DataVector& data2, const size_t linesToCompare, const double tolerance){
  int failure_found=EXIT_SUCCESS;

  if (data1.size() != data2.size() && linesToCompare > std::min(data1.size(),data2.size())){
    std::cout<< "two files have different numbers of lines needing comparing." <<std::endl;
    failure_found = VALUE_TOLERANCE_FAILURE_CODE;
  }

  if ( ! failure_found && linesToCompare > std::min(data1.size(),data2.size()) ) {
    std::cout<< "linesToCompare is greater than the number of data lines." <<std::endl;
    failure_found = VALUE_TOLERANCE_FAILURE_CODE;
  }

  if( ! failure_found )
  {
    for(size_t i=0; i< linesToCompare; ++i){
      std::vector< double > row1 = data1[i];
      std::vector< double > row2 = data2[i];
      const size_t column = row1.size();
      if (column != row2.size()){
        printf("In the line %zu, 2 data has different column.\n", i);
        failure_found = VALUE_TOLERANCE_FAILURE_CODE;
      }
      for (size_t j=0; j< column; ++j){
        if (std::abs(row1[j] - row2[j]) > tolerance){
          printf("In the line %zu, column %zu:\ndata1 has value: %.8f;\ndata2 has value: %.8f.\n", i+1,j+1,row1[j],row2[j]);
          printf("The difference is greater than the tolerance %g. \nSo data in two files are not close within specific tolerance.\nProgram exits.\n", tolerance);
          failure_found = VALUE_TOLERANCE_FAILURE_CODE;
        }
      }
    }
  }
  if ( ! failure_found )
  {
    printf("Two files are close within the specific tolerance %g \n", tolerance);
    printf("Program exits successfully.\n");
  }
  else{
    std::cout << "ERROR! ERROR! ERROR! The files do not match." << std::endl;
  }
  return failure_found;
}


int main(int argc, char * argv[])
{
  if( argc != 5 )
    {
    std::cerr << "USAGE: " << argv[0] << " <BFile1> <BFile2> <linesToCompare> <tolerance>" << std::endl;
    std::cerr << "Parameters have errors. Program exits." <<std::endl;
    return USAGE_FAILURE_CODE;
    }
  const int linesToCompare = atoi(argv[3]);
  const double tolerance = atof(argv[4]);


  //Open files to read data into vector
  char* buffer1 = NULL;
  if (0 != readFileIntoBuffer(argv[1], &buffer1)){
    return FILE_READ_FAILURE_CODE;
  }
  char* buffer2 = NULL;
  if (0 != readFileIntoBuffer(argv[2], &buffer2)){
    delete[] buffer1;
    return FILE_READ_FAILURE_CODE;
  }

  DataVector data1 = convertBufferIntoVector(buffer1);
  DataVector data2 = convertBufferIntoVector(buffer2);

  if (0 == compareRowColumn(data1, data2, linesToCompare, tolerance)){
    return EXIT_SUCCESS;
  }
  else{
    return EXIT_FAILURE;
  }
}
