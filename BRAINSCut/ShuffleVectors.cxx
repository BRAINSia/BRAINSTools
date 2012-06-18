#include "ShuffleVectors.h"
#include "BRAINSCutDataHandler.h"

#define MAX_LINE_SIZE 1000

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
    std::cout << "Error: Could not open ANN vector file"
              << std::endl
              << m_inputHeaderFilename;
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

  std::cout << "IVS = " << m_IVS << std::endl;
  std::cout << "OVS = " << m_OVS << std::endl;

  m_output_TVC = std::floor( (float)m_input_TVC * m_resampleProportion + 0.5F );
  std::cout << m_input_TVC << " * " << m_resampleProportion << " = "
            << "m_output_TVC == " << m_output_TVC << std::endl;
}

void
ShuffleVectors::WriteHeader()
{
  std::ofstream filestr;

  filestr.open( m_outputHeaderFilename.c_str() );
  if( !filestr.is_open() )
    {
    std::cout << "Error: Could not open ANN vector file"
              << std::endl;
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
ShuffleVectors::ShuffleOrder()
{
  typedef findUINT64Type<sizeof(unsigned long)>::unsigned64 unsigned64;

  vnl_random randgen;

#define randgen64()                                      \
  ( ( static_cast<unsigned64>( randgen.lrand32() ) << 32 ) \
    | static_cast<unsigned64>( randgen.lrand32() ) )

  std::ios::off_type *rval = new std::ios::off_type[m_output_TVC];

  if( rval == 0 )
    {
    std::cout << "Can't allocate shuffled ordering"
              << std::endl;
    }
  for( unsigned long i = 0; i < m_output_TVC; i++ )
    {
    rval[i] = static_cast<std::ios::off_type>( i );
    }
  // do the shuffle
  for( std::ios::off_type i = m_output_TVC - 1; i > 0; i-- )
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
  m_resampleProportion(0.0F)
{
}

ShuffleVectors::ShuffleVectors(const std::string& inputVectorFilename,
                               const std::string& outputVectorFilename,
                               float resampleProportion  ) :
  m_IVS(0),
  m_OVS(0),
  m_input_TVC(0),
  m_output_TVC(0),
  m_resampleProportion(0.0F)
{
  std::cout << "Shuffle Vectors of ======================================= " << std::endl
            << inputVectorFilename << " to " << std::endl
            << outputVectorFilename << std::endl
            << "========================================================== "
            << std::endl;

  if( inputVectorFilename == "" || outputVectorFilename == "" )
    {
    std::cout << " Filenames for inputVector and outputVector are neccessary"
              << std::endl;
    }
  else if( inputVectorFilename == outputVectorFilename )
    {
    std::cout << "ERROR:  Can not use the same file for input and output."
              << inputVectorFilename;
    }
  else
    {
    m_inputVectorFilename = inputVectorFilename;
    m_inputHeaderFilename = inputVectorFilename + ".hdr";
    m_outputVectorFilename = outputVectorFilename;
    m_outputHeaderFilename = outputVectorFilename + ".hdr";
    }

  if( !itksys::SystemTools::FileExists( inputVectorFilename.c_str(), false ) )
    {
    std::cout << "ERROR: Cannot open " << inputVectorFilename
              << ". \n The file does not exist."
              << std::endl;
    }
  m_resampleProportion = resampleProportion;
}

void
ShuffleVectors::Shuffling()
{
  // read binary
  std::ifstream inputVectorFileStream;

  inputVectorFileStream.open( m_inputVectorFilename.c_str(),
                              std::ios::in | std::ios::binary);
  if( !inputVectorFileStream.is_open() )
    {
    std::cout << "Can't open " << m_inputVectorFilename;
    }

  int shuffledFile =
    open(m_outputVectorFilename.c_str(),
         O_WRONLY | O_CREAT | O_TRUNC,
         S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH);
  if( shuffledFile == -1 )
    {
    std::cout << "Can't open output file " << m_outputVectorFilename;
    }
  // make a shuffled output ordering

  std::ios::off_type *randomOrder = ShuffleOrder();

  unsigned recordsize = ( m_IVS + m_OVS + 1 ) * sizeof( float );
  float *  buf = new float[m_IVS + m_OVS + 1];

  std::cout << "Writing a shuffled output file "
            << m_outputVectorFilename
            << std::endl;

  float VectorsPerPercent = m_input_TVC / 100.0F;
  int   current_percent = 0;
  /**
   * shuffle process goes through entire input vector file and then
   * include only those have (randomly assigned) index less
   * than a desired output size
   */
  for( unsigned int vectorIndex = 0;
       vectorIndex < m_output_TVC;
       vectorIndex++ )
    {
    if( vectorIndex > current_percent * VectorsPerPercent )
      {
      current_percent += 1;
      if( current_percent % 5 == 0 )
        {
        std::cout << current_percent << "% " << std::endl;
        }
      }
    if( inputVectorFileStream.eof() )
      {
      if( (vectorIndex + 1) % m_input_TVC == 1 ) // vector index starts from 0
        {
        std::cout << "*** Re-open the vector file" << std::endl;
        // read input vector file stream from the first again
        inputVectorFileStream.close();
        inputVectorFileStream.open( m_inputVectorFilename.c_str(),
                                    std::ios::in | std::ios::binary);
        }
      else
        {
        std::cerr << "Premature end of file at record "
                  << vectorIndex << std::endl
                  << (vectorIndex + 1) % m_input_TVC << " != 1" << std::endl;
        break; // TODO throw error here
        }
      }
    // read in the record
    inputVectorFileStream.read( (char *)buf, recordsize );
    /*
    std::cout<< " randomOrder[ "<< vectorIndex<<" ] :: "
             << randomOrder[vectorIndex]
             << " < " << m_output_TVC
             << " ? " << std::endl;
             */
    if( randomOrder[vectorIndex] < static_cast<std::ios::off_type>( m_output_TVC) )
      {
      std::ios::off_type seekval =
        randomOrder[vectorIndex] * static_cast<std::ios::off_type>( recordsize );
      lseek(shuffledFile, seekval, SEEK_SET);
      (void)write(shuffledFile, (const char *) buf, recordsize);
      /*
      // debugging code
      for( int dummy_i = 0; dummy_i< m_IVS  + m_OVS +1; dummy_i++ )
        {
        std::cout<<buf[dummy_i]<<" ";
        }
      std::cout<< " randomOrder[ "<< vectorIndex<<" ] :: "
               << randomOrder[vectorIndex]
               << " @ "
               <<seekval <<std::endl;
      //*/
      }

    if( buf[m_IVS  + m_OVS] != LineGuard )
      {
      std::cerr << "Record not properly terminated by sentinel value ::  "
                << buf[m_IVS  + m_OVS] << " != "
                << LineGuard
                << " at Vector index " << vectorIndex
                << std::endl;
      exit(EXIT_FAILURE); // TODO throw error here
      }
    }
  close(shuffledFile);
  std::cout << "done." << std::endl;
  delete[] buf;
  delete[] randomOrder;
}
