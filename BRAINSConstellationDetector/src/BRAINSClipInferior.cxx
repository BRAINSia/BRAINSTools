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
/*
 * Author: Hans J. Johnson
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * Redistributions of source code must retain the above copyright notice, this
 * list of conditions and the following disclaimer.
 *
 * Redistributions in binary form must reproduce the above copyright notice,
 * this list of conditions and the following disclaimer in the documentation
 * and/or other materials provided with the distribution.
 *
 * Neither the name of the Nathan Kline Institute nor the names of its
 * contributors may be used to endorse or promote products derived from this
 * software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 */
#include "BRAINSThreadControl.h"
#include "BRAINSClipInferiorCLP.h"
#include "itkIO.h"
#include "landmarksConstellationCommon.h"

#include "ChopImageBelowLowerBound.h"

//
//
// ////////////////////////////////////////////////////////////////////////////////////////////////
int
main(int argc, char * argv[])
{
  PARSE_ARGS;
  BRAINSRegisterAlternateIO();
  const BRAINSUtils::StackPushITKDefaultNumberOfThreads TempDefaultNumberOfThreadsHolder(numberOfThreads);
  bool                                                  verbose = true;

  std::cout << "================================================================" << std::endl;
  std::cout << "Processing: " << inputVolume << std::endl;

  // /////////////////////////////////////////////////////////////////////////////////////////////
  if (verbose)
  {
    printf("-------------------------------------------------------\n");
    printf("inputVolume: %s\n", inputVolume.c_str());
    printf("outputVolume: %s\n", outputVolume.c_str());
    printf("acLowerBound: %f\n", acLowerBound);
  }
  if (outputVolume.empty())
  {
    std::cout << "ERROR:  Missing output file name." << std::endl;
    std::cout << "        Please specify -o <filename> or --outputVolume <filename>" << std::endl;
    return -1;
  }
  // /////////////////////////////////////////////////////////////////////////////////////////////
  short BackgroundFillValue;
  if (backgroundFillValueString == std::string("BIGNEG"))
  {
    BackgroundFillValue = -32768;
  }
  else
  {
    BackgroundFillValue = std::stoi(backgroundFillValueString.c_str());
  }
  // //////////////////////////////////////////////////////////////////////////
  SImageType::Pointer image = itkUtil::ReadImage<SImageType>(inputVolume);
  if (image.IsNull())
  {
    printf("\nCould not open image %s, aborting ...\n\n", inputVolume.c_str());
    return EXIT_FAILURE;
  }

  // we need a DOUBLE constant, not a FLOAT constant, for exact switch
  // comparisons.
  constexpr double thousand = 1000.0;
  if (acLowerBound != thousand)
  {
    double PhysicalLowerBound = /* ACy when zero-centered is ... */ 0.0 - acLowerBound;
    ChopImageBelowLowerBound<SImageType>(image, BackgroundFillValue, PhysicalLowerBound);
  }
  itkUtil::WriteImage<SImageType>(image, outputVolume);
  if (verbose)
  {
    printf("Wrote outputVolume: %s\n", outputVolume.c_str());
    printf("-------------------------------------------------------\n");
  }
  return 0;
}
