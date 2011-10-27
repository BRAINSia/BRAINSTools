#include "ShuffleVectors.h"
#include "ShuffleVectorsModuleCLP.h"

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Shuffling Vectors
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
int
main(int argc, char * *argv)
{
  PARSE_ARGS;
  if( argc < 4 )
    {
    std::cerr << "USAGE::::" << std::endl
              << argv[0] << std::endl
              << " [inputVectorFileBaseName] [outputVectorFileBaseName] [downsamplesize] "
              << std::endl;
    std::cerr
      <<
      "downsample size of 1 will be the same size as the input images, downsample size of 3 will throw 2/3 the vectors away."
      << std::endl;
    exit(1);
    }
  // Shuffled the vector:
  ShuffleVectors * my_ShuffleVector = new ShuffleVectors(  std::string( argv[1] ),
                                                           std::string( argv[2] ),
                                                           atoi( argv[3] ) );
  my_ShuffleVector->ReadHeader();
  my_ShuffleVector->Shuffling();
  my_ShuffleVector->WriteHeader();

  return 0;
}
