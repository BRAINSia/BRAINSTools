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
#include "ShuffleVectors.h"
#include "ShuffleVectorsModuleCLP.h"
#include <BRAINSCommonLib.h>

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Shuffling Vectors
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
int
main(int argc, char ** argv)
{
  PARSE_ARGS;
  BRAINSRegisterAlternateIO();
  // Shuffled the vector:
  ShuffleVectors * my_ShuffleVector =
    new ShuffleVectors(inputVectorFileBaseName, outputVectorFileBaseName, resampleProportion);
  my_ShuffleVector->ReadHeader();
  my_ShuffleVector->Shuffling();
  my_ShuffleVector->WriteHeader();

  return EXIT_SUCCESS;
}
