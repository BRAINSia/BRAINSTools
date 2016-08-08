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
#include <BRAINSCommonLib.h>
#include <itkFFTWCommon.h>

void BRAINSRegisterAlternateIO(void)
{
}

//This is intended to be called one time
void FFTWInit(const std::string path_for_wisdom)
{
  //Environmental variables
  //itksys::SystemTools::GetEnv("ITK_FFTW_PLAN_RIGOR", "STRING");
  //ITK_FFTW_PLAN_RIGOR   - Defines how aggressive the generation of wisdom should be.
  //         FFTW_ESTIMATE - super fast guess to plan, and often mediocre performance for odd sizes
  //         FFTW_MEASURE - quick planning wiht heuristics, and usually good performance
  //         FFTW_PATIENT -    slow planning with heuristics, and almost always best performance (skip testing
  // odd-ball cases)
  //         FFTW_EXHAUSTIVE - very slow planning by checking odd-ball cases, often does not produce better results
  //ITK_FFTW_READ_WISDOM_CACHE  - Defines if a wisdom file cache should
  //                              be read if found.  (it is "On" by default)
  //ITK_FFTW_WRITE_WISDOM_CACHE  - Defines if generated wisdom file cache
  //                               should be written (it is "Off" by default)
  //ITK_FFTW_WISDOM_CACHE_BASE - Defines the base directory where the
  //                             fftw wisdom cache will be placed,
  //                             this is intended to be used with auto-
  //                             generated cache file names
  //ITK_FFTW_WISDOM_CACHE_FILE - Defines the full name of the cache
  //                             file to be generated.  If this is
  //                             set, then ITK_FFTW_WISDOM_CACHE_BASE
  //                             is ignored.
  if( path_for_wisdom.length() > 0) // If empty, just use the default.
  {
    itk::FFTWGlobalConfiguration::SetWisdomCacheBase(path_for_wisdom);
  }
  itk::FFTWGlobalConfiguration::SetReadWisdomCache(true);
  itk::FFTWGlobalConfiguration::SetWriteWisdomCache(true);
  itk::FFTWGlobalConfiguration::ImportDefaultWisdomFileDouble();
  itk::FFTWGlobalConfiguration::ImportDefaultWisdomFileFloat();
  itk::FFTWGlobalConfiguration::SetPlanRigor(FFTW_EXHAUSTIVE);
  //itk::FFTWGlobalConfiguration::SetPlanRigor(FFTW_MEASURE);
  //std::string temp = itk::FFTWGlobalConfiguration::GetWisdomFileDefaultBaseName();
  //FFTW_MEASURE, FFTW_PATIENT, or FFTW_EXHAUSTIVE
  // We need to ensure that we read in the default wisdom files from the default locations
  itk::FFTWGlobalConfiguration::ImportDefaultWisdomFileFloat();
  itk::FFTWGlobalConfiguration::ImportDefaultWisdomFileDouble();
  std::cout << itk::FFTWGlobalConfiguration::GetWisdomCacheBase() << std::endl;
  std::cout << itk::FFTWGlobalConfiguration::GetWisdomFileDefaultBaseName() << std::endl;
}
