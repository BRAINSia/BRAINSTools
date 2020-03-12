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

#include "itkIO.h"
#include "BRAINSTrimForegroundInDirectionCLP.h"
#include "TrimForegroundInDirection.h"
#include "landmarkIO.h"
#include "BRAINSThreadControl.h"

//
//
// ////////////////////////////////////////////////////////////////////////////////////////////////
int
main(int argc, char * argv[])
{
  // file pointer for opening the setup file
  // /////////////////////////////////////////////////////////////////////////////////////////////
  PARSE_ARGS;
  BRAINSRegisterAlternateIO();
  const BRAINSUtils::StackPushITKDefaultNumberOfThreads TempDefaultNumberOfThreadsHolder(numberOfThreads);
  std::cout << "================================================================" << std::endl;
  std::cout << "Processing: " << inputVolume << std::endl;

  // //////////////////////////////////////////////////////////////////////////
  SImageType::Pointer volOrig = itkUtil::ReadImage<SImageType>(inputVolume);
  if (volOrig.IsNull())
  {
    std::cerr << "Could not open image " << inputVolume.c_str() << std::endl;
    return EXIT_FAILURE;
  }

  // /////////////////////////////////////////////////////////////////////////////////////////////
  if (directionCode == 0)
  {
    std::cout << "Your directionCode of 0 has selected the program default of 3 (maximize Superior/Inferior)."
              << std::endl;
    directionCode = 3;
  }
  unsigned int axis = itk::Math::abs(directionCode) - 1;
  if (axis > 2)
  {
    std::cout << "Your directionCode was too large so we will use the program default, axis 2 (Superior/Inferior)."
              << std::endl;
    axis = 2;
  }

  // /////////////////////////////////////////////////////////////////////////////////////////////
  SImageType::PixelType BackgroundFillValue = 0;
  if (backgroundFillValueString == std::string("BIGNEG"))
  {
    BackgroundFillValue = -32768;
  }
  else
  {
    BackgroundFillValue = std::stoi(backgroundFillValueString.c_str());
  }

  // /////////////////////////////////////////////////////////////////////////////////////////////
  SImageType::Pointer volResult;

  TrimForegroundInDirection(
    volResult, volOrig, axis, otsuPercentileThreshold, closingSize, headSizeLimit, BackgroundFillValue);
  itkUtil::WriteImage<SImageType>(volResult, outputVolume);
  return 0;
}
