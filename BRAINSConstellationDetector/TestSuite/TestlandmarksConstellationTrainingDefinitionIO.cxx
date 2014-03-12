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
#include "../src/landmarksConstellationTrainingDefinitionIO.h"

int main(int argc, char * *argv)
{
  if( argc != 2 )
    {
    std::cerr << "Usage: testAcpcmodelSetupClass <model txt file>"
              << std::endl;
    std::cerr.flush();
    return EXIT_FAILURE;
    }
  std::string                                filename(argv[1]);
  landmarksConstellationTrainingDefinitionIO m;
  if( m.ReadFile(filename) == -1 )
    {
    std::cerr << "Error reading " << filename
              << std::endl;
    std::cerr.flush();
    return EXIT_FAILURE;
    }
  std::cout << m;
  return EXIT_SUCCESS;
}
