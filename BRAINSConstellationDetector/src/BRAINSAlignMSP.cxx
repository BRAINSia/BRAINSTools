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
 * SUBidentityTransformITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, identityTransformRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 */
#include <sstream>
#include <cmath>
#include <itkIntensityWindowingImageFilter.h>
#include "BRAINSThreadControl.h"
#include "landmarksConstellationCommon.h"
#include "itkFindCenterOfBrainFilter.h"
#include "itkIO.h"
#include "BRAINSAlignMSPCLP.h"
#include "GenericTransformImage.h"

int main(int argc, char *argv[])
{
  std::cout.precision(10);

  // /////////////////////////////////////////////////////////////////////////////////////////////
  PARSE_ARGS;
  BRAINSRegisterAlternateIO();

  const BRAINSUtils::StackPushITKDefaultNumberOfThreads TempDefaultNumberOfThreadsHolder(numberOfThreads);

  LMC::globalverboseFlag = verbose;

  globalResultsDir = resultsDir;
  globalImagedebugLevel = writedebuggingImagesLevel;
  // /////////////////////////////////////////////////////////////////////////////////////////////
  // read information from the setup file, allocate some memories, and
  // initialize some variables

  // Since these are oriented images, the reorientation should not be necessary.
  SImageType::Pointer volOrig = itkUtil::ReadImage<SImageType>(inputVolume);
  if( volOrig.IsNull() )
    {
    printf( "\nCould not open image %s, aborting ...\n\n", inputVolume.c_str() );
    return EXIT_FAILURE;
    }
  SImageType::Pointer image;

  if( rescaleIntensities == true )
    {
    itk::StatisticsImageFilter<SImageType>::Pointer stats =
      itk::StatisticsImageFilter<SImageType>::New();
    stats->SetInput(volOrig);
    stats->Update();
    SImageType::PixelType minPixel( stats->GetMinimum() );
    SImageType::PixelType maxPixel( stats->GetMaximum() );

    if( trimRescaledIntensities > 0.0 )
      {
      // REFACTOR: a histogram would be traditional here, but seems
      // over-the-top;
      // I did this because it seemed to me if I knew mean, sigma, max and min,
      // then I know Something about extreme outliers.

      double meanOrig( stats->GetMean() );
      double sigmaOrig( stats->GetSigma() );

      // REFACTOR:  In percentiles, 0.0005 two-tailed has worked in the past.
      // It only makes sense to trim the upper bound since the lower bound would
      // most likely
      // represent a large region of air around the head.  But this is not so
      // when using a mask.
      // For one-tailed, an error of 0.001 corresponds to 3.29052 standard
      // deviations of normal.
      // For one-tailed, an error of 0.0001 corresponds to 3.8906 standard
      // deviations of normal.
      // For one-tailed, an error of 0.00001 corresponds to 4.4172 standard
      // deviations of normal.
      // Naturally, the constant should default at the command line, ...

      double variationBound( ( maxPixel - meanOrig ) / sigmaOrig );
      double trimBound(variationBound - trimRescaledIntensities);
      if( trimBound > 0.0 )
        {
        maxPixel = static_cast<SImageType::PixelType>( maxPixel - trimBound * sigmaOrig );
        }
      }

    itk::IntensityWindowingImageFilter<SImageType, SImageType>::Pointer remapIntensityFilter =
      itk::IntensityWindowingImageFilter<SImageType, SImageType>::New();
    remapIntensityFilter->SetInput(volOrig);
    remapIntensityFilter->SetOutputMaximum(rescaleIntensitiesOutputRange[1]);
    remapIntensityFilter->SetOutputMinimum(rescaleIntensitiesOutputRange[0]);
    remapIntensityFilter->SetWindowMinimum(minPixel);
    remapIntensityFilter->SetWindowMaximum(maxPixel);
    remapIntensityFilter->Update();

    image = remapIntensityFilter->GetOutput();
    }
  else
    {
    image = volOrig;
    }

  // Find center of head mass
  std::cout << "\nFinding center of head mass..." << std::endl;
  typedef itk::FindCenterOfBrainFilter<SImageType>                        FindCenterFilter;
  FindCenterFilter::Pointer findCenterFilter = FindCenterFilter::New();
  findCenterFilter->SetInput(image);
  findCenterFilter->SetAxis(2);
  findCenterFilter->SetOtsuPercentileThreshold(0.01);
  findCenterFilter->SetClosingSize(7);
  findCenterFilter->SetHeadSizeLimit(700);
  findCenterFilter->SetBackgroundValue(0);
  findCenterFilter->Update();
  SImagePointType centerOfHeadMass = findCenterFilter->GetCenterOfBrain();

  RigidTransformType::Pointer Tmsp = RigidTransformType::New();
  ComputeMSP_Easy(image, Tmsp, centerOfHeadMass, mspQualityLevel);

  // /////////////////////////////////////////////////////////////////////////////////////////////
  short BackgroundFillValue;
  if( backgroundFillValueString == std::string("BIGNEG") )
    {
    BackgroundFillValue = -32768;
    }
  else
    {
    BackgroundFillValue = atoi( backgroundFillValueString.c_str() );
    }

    {
    // Remember:  the Data is Moving's, the shape is Fixed's.
    SImageType::Pointer interpImage = TransformResample<SImageType, SImageType>(
        image.GetPointer(), image.GetPointer(), BackgroundFillValue,
        GetInterpolatorFromString<SImageType>(interpolationMode).GetPointer(), Tmsp.GetPointer() );
    itkUtil::WriteImage<SImageType>(interpImage, resampleMSP);
    }
  if( globalImagedebugLevel > 3 )
    {
    const std::string ORIG_ImagePlane( globalResultsDir + "/ORIG_PLANE_" + itksys::SystemTools::GetFilenameName(
                                         inputVolume) );
    CreatedebugPlaneImage(image, Tmsp, ORIG_ImagePlane);
    }
  // TODO:  Add more features to this program:
  //       1) Resample to a given space (256^3, 1.0mm^3)
  //       2) Output using Windowed Sinc for best results.
  return 0;
}
