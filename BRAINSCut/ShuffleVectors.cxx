#include "ShuffleVectors.h"

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Usage Example
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
/*
main(int argc, char **argv)
{
  if ( argc < 4 )
  {
    std::cerr << "USAGE::::" <<std::endl
              << argv[0] <<std::endl
              << " [inputVectorFileBaseName] [outputVectorFileBaseName] [downsamplesize] "
              << std::endl;
    std::cerr
    <<
    "downsample size of 1 will be the same size as the input images, downsample size of 3 will throw 2/3 the vectors away."
    << std::endl;

  }
  //Shuffled the vector:
  ShuffleVectors * my_ShuffleVector = new ShuffleVectors (  std::string( argv[1] ),
                                                          std::string( argv[2] ),
                                                          atoi( argv[3] ) );
  my_ShuffleVector -> ReadHeader();
  my_ShuffleVector -> Shuffling();
  my_ShuffleVector -> WriteHeader();


  return 0;
}
*/

void
ShuffleVectors::ReadHeader()
{
  std::ifstream filestr;

  filestr.open( m_inputHeaderFilename.c_str() );
  if( !filestr.is_open() )
    {
    itkGenericExceptionMacro(<< "Error: Could not open ANN vector file"
                             << std::endl
                             << m_inputHeaderFilename);
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
      iss >> m_IVS;
      }
    else if( temp == "OVS" )
      {
      iss >> m_OVS;
      }
    else if( temp == "TVC" )
      {
      iss >> m_input_TVC;
      }
    }

  filestr.close();

  m_output_TVC = m_input_TVC / m_downSampleSize;
}

void
ShuffleVectors::WriteHeader()
{
  std::ofstream filestr;

  filestr.open( m_outputHeaderFilename.c_str() );
  if( !filestr.is_open() )
    {
    itkGenericExceptionMacro(<< "Error: Could not open ANN vector file")
    }
  filestr << "IVS " <<  m_IVS << std::endl;
  filestr << "OVS " <<  m_OVS << std::endl;
  filestr << "TVC " <<  m_output_TVC << std::endl;
  filestr << "SHUFFLED TRUE" << std::endl;
  filestr.close();
}

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

std::ios::off_type *
ShuffleVectors::ShuffledOrder()
{
  typedef findUINT64Type<sizeof(unsigned long)>::unsigned64 unsigned64;

  vnl_random randgen;

#define randgen64()                                      \
  ( ( static_cast<unsigned64>( randgen.lrand32() ) << 32 ) \
    | static_cast<unsigned64>( randgen.lrand32() ) )

  std::ios::off_type *rval = new std::ios::off_type[m_input_TVC];

  if( rval == 0 )
    {
    itkGenericExceptionMacro(<< "Can't allocate shuffled ordering")
    }
  for( unsigned long i = 0; i < m_input_TVC; i++ )
    {
    rval[i] = static_cast<std::ios::off_type>( i );
    }
  // do the shuffle
  for( std::ios::off_type i = m_input_TVC - 1; i > 0; i-- )
    {
    std::ios::off_type j( randgen64() % ( i + 1 ) ); // rand() % (i+1);
    std::ios::off_type tmp(rval[i]);
    rval[i] = rval[j];
    rval[j] = tmp;
    }
  return rval;
}

//
// Constructors
//
ShuffleVectors::ShuffleVectors() :
  m_IVS(0),
  m_OVS(0),
  m_input_TVC(0),
  m_output_TVC(0),
  m_downSampleSize(0)
{
}

ShuffleVectors::ShuffleVectors(const std::string& inputVectorFilename, const std::string& outputVectorFilename,
                               int downSampleSize  ) :
  m_IVS(0),
  m_OVS(0),
  m_input_TVC(0),
  m_output_TVC(0),
  m_downSampleSize(0)
{
  std::cout << "Shuffle Vectors of ======================================= " << std::endl
            << inputVectorFilename << " to " << std::endl
            << outputVectorFilename << std::endl
            << "========================================================== "
            << std::endl;

  if( inputVectorFilename == "" || outputVectorFilename == "" )
    {
    itkGenericExceptionMacro(<< " Filenames for inputVector and outputVector are neccessary");
    }
  else if( inputVectorFilename == outputVectorFilename )
    {
    itkGenericExceptionMacro(<< "ERROR:  Can not use the same file for input and output."
                             << inputVectorFilename);
    }
  else
    {
    m_inputVectorFilename = inputVectorFilename;
    m_inputHeaderFilename = inputVectorFilename + ".hdr";
    m_outputVectorFilename = outputVectorFilename;
    m_outputHeaderFilename = outputVectorFilename + ".hdr";
    }

  m_downSampleSize = downSampleSize;
}

void
ShuffleVectors::Shuffling()
{
  // read binary
  std::ifstream binfile;

  binfile.open( m_inputVectorFilename.c_str(),
                std::ios::in | std::ios::binary);
  if( !binfile.is_open() )
    {
    itkGenericExceptionMacro(<< "Can't open " << m_inputVectorFilename );
    }

  int shuffledFile =
    open(m_outputVectorFilename.c_str(),
         O_WRONLY | O_CREAT | O_TRUNC,
         S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH);
  if( shuffledFile == -1 )
    {
    itkGenericExceptionMacro(<< "Can't open output file " << m_outputVectorFilename);
    }
  // make a shuffled output ordering
  std::cout.flush();

  std::ios::off_type *order = ShuffledOrder();

  unsigned recordsize = ( m_IVS + m_OVS + 1 ) * sizeof( float );
  float *  buf = new float[m_IVS + m_OVS + 1];

  std::cout << "Writing a shuffled output file "
            << m_outputVectorFilename
            << std::endl;
  std::cout.flush();

  float VectorsPerPercent = m_output_TVC / 100.0F;
  int   current_percent = 0;
  for( unsigned int i = 0; i < m_output_TVC; i++ )
    {
    if( i > current_percent * VectorsPerPercent )
      {
      current_percent += 1;
      if( current_percent % 5 == 0 )
        {
        std::cout << current_percent << "% ";
        std::cout.flush();
        }
      }
    if( binfile.eof() )
      {
      std::cerr << "Premature end of file at record "
                << i << std::endl;
      break; // TODO throw error here
      }
    // read in the record
    binfile.read( (char *)buf, recordsize );

    if( buf[m_IVS  + m_OVS] != AUTOSEG_VEC_SENTINEL )
      {
      std::cerr << "Record not properly terminated by sentinel value ::  "
                << buf[m_IVS  + m_OVS]
                << std::endl;
      break; // TODO throw error here
      }
    // Now only write out vector if it is supposed to be placed in the smaller
    // image size.
    if( order[i] < static_cast<std::ios::off_type>( m_output_TVC) )
      {
      std::ios::off_type seekval =
        order[i] * static_cast<std::ios::off_type>( recordsize );
      lseek(shuffledFile, seekval, SEEK_SET);
      (void)write(shuffledFile, buf, recordsize);
      }
    }
  close(shuffledFile);
  std::cout << "done." << std::endl;
  delete[] buf;
  delete[] order;
}
