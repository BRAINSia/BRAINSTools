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
#include <stdint.h>

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
  ShuffleVectors( const std::string& inputFilename, const std::string& outputFilename, int downSampleSize = 1 );
  ~ShuffleVectors();

  void Shuffling();

  void ReadHeader();

  void WriteHeader();

private:
  std::string TempName( const char *s );

  std::ios::off_type * ShuffleOrder();

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
  unsigned int m_downSampleSize;
};

#endif
