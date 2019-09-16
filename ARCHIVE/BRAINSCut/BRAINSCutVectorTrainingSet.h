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
#ifndef BRAINSCutVectorTrainingSet_h
#define BRAINSCutVectorTrainingSet_h

#include "BRAINSCutDataHandler.h"

// using neuralNetType = CvANN_MLP_Revision;
using neuralNetType = cv::ml::ANN_MLP;

namespace
{
template <size_t LongSize>
class findUINT64Type
{};
template <>
class findUINT64Type<4>
{
public:
  using unsigned64 = unsigned long long;
};
template <>
class findUINT64Type<8>
{
public:
  using unsigned64 = unsigned long;
};
} // namespace

class BRAINSCutVectorTrainingSet
{
public:
  /* constructor
   * :: the paried training set should be constructed with vector file name.
   */
  BRAINSCutVectorTrainingSet(std::string vectorFilenamePrefix);
  ~BRAINSCutVectorTrainingSet();

  static const unsigned int MAXIMUMCHAR = 100;

  void
  ReadHeaderFileInformation();

  void
  SetRecordSize(); // TODO : name is subject change. Keep it for now consistency to the previous version

  unsigned int
  GetRecordSize();

  void
  SetBufferRecordSize();

  void
  SetShuffled(bool shuffled);

  int
  GetTotalVectorSize();

  int
  GetInputVectorSize();

  int
  GetOutputVectorSize();

  void
  PrintDebuggingMessage(std::string msg);

  scalarType *
  ReadBufferFromFileStream(std::ifstream & fileStream);

  void
  RandomizeTrainingVector();

  pairedTrainingSetType *
  GetTrainingDataSet();

  void
  SetTrainingSubSet(unsigned int count);

  pairedTrainingSetType *
  GetTrainingSubSet(unsigned int count);

  pairedTrainingSetType *
  DownSampleTrainingDataSet(const unsigned int subSampleSize);

  void
  WriteVectorFile();

  void
  SetNumberOfSubSet(const unsigned int count = 1);

  // INFO: REGINA  All "Get" functions should be const
  unsigned int
  GetNumberOfSubSet();

  // TODO
  // IMPLEMEMNT THIS::
  void
  SaveCurrentSubSet(std::string filename);

private:
  // INFO: REGINA these all need to have "m_" prefix
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

  unsigned int currentSubSetID; // goes from 0,1,..
};
#endif
