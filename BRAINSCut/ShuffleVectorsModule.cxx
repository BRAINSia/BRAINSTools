#include "ShuffleVectors.h"
#include "ShuffleVectorsModuleCLP.h"

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Shuffling Vectors
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
int
main(int argc, char * *argv)
{
  PARSE_ARGS;
  // Shuffled the vector:
  ShuffleVectors * my_ShuffleVector = new ShuffleVectors(  inputVectorFileBaseName,
                                                           outputVectorFileBaseName,
                                                           resampleProportion);
  my_ShuffleVector->ReadHeader();
  my_ShuffleVector->Shuffling();
  my_ShuffleVector->WriteHeader();

  return EXIT_SUCCESS;
}
