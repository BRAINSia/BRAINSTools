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
#ifndef ShuffleVectors_h
#define ShuffleVectors_h

#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vnl/vnl_random.h>
#include <itksys/SystemTools.hxx>

class ShuffleVectors
{
  //
  // Description::
  // ShuffleVectors Class modified from ShuffleVector Function.
  // The program will shuffle the 'inputFilename' Vector and store as
  //    'outputFilename' vector with proper header file.
  //    Only assumption is that header file name s to be
  // 'AnyVectorFilename'+'.hdr'.
  //    The output header file will be flaged as shuffle = True
  //    - Eun Young (Regina) Kim
public:
  ShuffleVectors();
  ShuffleVectors( const std::string& inputFilename, const std::string& outputFilename, float downSampleSize = 1.0F );
  ~ShuffleVectors();

  void Shuffling();

  void ReadHeader();

  void WriteHeader();

private:
  std::string TempName( const char *s );

  void ShuffleOrder(std::vector<std::ios::off_type> & rval) const;

  //
  // Member Variables::
  //
  // - Filename ::
  std::string m_inputVectorFilename;
  std::string m_inputHeaderFilename;
  std::string m_outputVectorFilename;
  std::string m_outputHeaderFilename;
  //
  // - Header Info ::
  //
  unsigned int  m_IVS;
  unsigned int  m_OVS;
  unsigned long m_input_TVC;
  unsigned long m_output_TVC;
  //
  // - Down Sampling Size
  //
  float m_resampleProportion;
};

#endif
