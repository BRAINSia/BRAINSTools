#ifndef BRAINSCutVectorTrainingSet_h
#define BRAINSCutVectorTrainingSet_h

#include "BRAINSCutPrimary.h"

typedef CvANN_MLP_Revision neuralNetType;

namespace
{
template <size_t LongSize>
class findUINT64Type
{
};
template <>
class findUINT64Type<4>
{
public: typedef unsigned long long unsigned64;
};
template <>
class findUINT64Type<8>
{
public:  typedef unsigned long unsigned64;
};
}

class BRAINSCutVectorTrainingSet
{
public:
  /* constructor
   * :: the paried training set should be constructed with vector file name.
   */
  BRAINSCutVectorTrainingSet(std::string vectorFilenamePrefix);

  static const unsigned int MAXIMUMCHAR = 100;

  void         ReadHeaderFileInformation();

  void         SetRecordSize();   // TODO : name is subject change. Keep it for now consistency to the previous version

  unsigned int GetRecordSize();

  void SetBufferRecordSize();

  int GetTotalVectorSize();

  int GetInputVectorSize();

  int GetOutputVectorSize();

  void                    PrintDebuggingMessage(std::string msg);

  scalarType *            ReadBufferFromFileStream( std::ifstream& fileStream );

  std::ios::off_type *    ShufflingOrder();

  void                    ShuffleVectors();

  pairedTrainingSetType * GetTrainingDataSet();

  void                    SetTrainingSubSet( unsigned int count );

  pairedTrainingSetType * GetTrainingSubSet( unsigned int count );

  pairedTrainingSetType * DownSampleTrainingDataSet( const unsigned int subSampleSize );

  void                    WriteVectorFile();

  void                    SetNumberOfSubSet( const unsigned int count = 1 );

  unsigned int            GetNumberOfSubSet();

  // TODO
  // IMPLEMEMNT THIS::
  void                    SaveCurrentSubSet( std::string filename );

private:
  /** file names */
  std::string trainingVectorFilename;
  std::string trainingHeaderFilename;

  /** size information from Header*/
  int totalVectorSize;
  int inputVectorSize;
  int outputVectorSize;

  bool               shuffled;
  std::ios::off_type recordSize;
  unsigned int       bufferRecordSize;

  /* input/output vector */
  unsigned int            numberOfSubSet;
  pairedTrainingSetType * currentTrainingSubSet;

  unsigned int currentSubSetID;   // goes from 0,1,..
};
#endif
