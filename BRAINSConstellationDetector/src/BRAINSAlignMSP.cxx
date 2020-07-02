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

int
main(int argc, char * argv[])
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
  using LOCALImageType = SImageType;
  using FloatImageType = itk::Image<float, 3>;
  FloatImageType::Pointer volOrig = itkUtil::ReadImage<FloatImageType>(inputVolume);
  if (volOrig.IsNull())
  {
    printf("\nCould not open image %s, aborting ...\n\n", inputVolume.c_str());
    return EXIT_FAILURE;
  }

  LOCALImageType::Pointer internal_image = [=]() -> LOCALImageType::Pointer {
    itk::RescaleIntensityImageFilter<FloatImageType, LOCALImageType>::Pointer remapIntensityFilter =
      itk::RescaleIntensityImageFilter<FloatImageType, LOCALImageType>::New();
    remapIntensityFilter->SetInput(volOrig);
    remapIntensityFilter->SetOutputMaximum(std::numeric_limits<LOCALImageType::PixelType>::max());
    remapIntensityFilter->SetOutputMinimum(std::numeric_limits<LOCALImageType::PixelType>::min());
    remapIntensityFilter->Update();
    return remapIntensityFilter->GetOutput();
  }();

  RigidTransformType::Pointer orig2msp_img_tfm = nullptr;

  SImagePointType orig_lmk_CenterOfHeadMass;
  // load corresponding landmarks from file, and define MSP from those landmarks.
  const LandmarksMapType constant_orig_lmks =
    (!LandmarkPoints.empty()) ? ReadSlicer3toITKLmk(LandmarkPoints) : LandmarksMapType{};
  if ((constant_orig_lmks.find("AC") != constant_orig_lmks.end()) &&
      (constant_orig_lmks.find("PC") != constant_orig_lmks.end()) &&
      (constant_orig_lmks.find("RP") != constant_orig_lmks.end()))
  {
    orig2msp_img_tfm = GetACPCAlignedZeroCenteredTransform(constant_orig_lmks);
  }
  else
  {
    // Estimate the center of head mass
    std::cout << "\nFinding center of head mass..." << std::endl;
    using FindCenterFilter = itk::FindCenterOfBrainFilter<LOCALImageType>;
    FindCenterFilter::Pointer findCenterFilter = FindCenterFilter::New();
    findCenterFilter->SetInput(internal_image);
    findCenterFilter->SetAxis(2);
    findCenterFilter->SetOtsuPercentileThreshold(0.01);
    findCenterFilter->SetClosingSize(7);
    findCenterFilter->SetHeadSizeLimit(700);
    findCenterFilter->SetBackgroundValue(0);
    findCenterFilter->Update();
    orig_lmk_CenterOfHeadMass = findCenterFilter->GetCenterOfBrain();
    double cc = -123.456;
    orig2msp_img_tfm = ComputeMSP(internal_image, orig_lmk_CenterOfHeadMass, mspQualityLevel, cc);
  }

  if (!resampleMSPLandmarkPoints.empty() && !LandmarkPoints.empty())
  {
    RigidTransformType::Pointer orig2msp_lmk_tfm = RigidTransformType::New();
    orig2msp_img_tfm->GetInverse(orig2msp_lmk_tfm);
    if (orig2msp_lmk_tfm.IsNull())
    {
      std::cerr << "The input transform not invertable." << std::endl;
      return EXIT_FAILURE;
    }

    LandmarksMapType transformedLandmarks;
    for ( const auto & elem : constant_orig_lmks)
    {
      transformedLandmarks[elem.first] = orig2msp_lmk_tfm->TransformPoint(elem.second);
    }

    WriteITKtoSlicer3Lmk(resampleMSPLandmarkPoints, transformedLandmarks);
    std::cout << "The transformed landmarks file is written." << std::endl;
  }

  // /////////////////////////////////////////////////////////////////////////////////////////////
  const short BackgroundFillValue = [=]() -> short {
    if (backgroundFillValueString == std::string("BIGNEG"))
    {
      return -32768;
    }
    return std::stoi(backgroundFillValueString.c_str());
  }();

  //Now write out the resampled image
  {
    LOCALImageType::Pointer image;
    if (rescaleIntensities)
    {
      itk::StatisticsImageFilter<FloatImageType>::Pointer stats = itk::StatisticsImageFilter<FloatImageType>::New();
      stats->SetInput(volOrig);
      stats->Update();
      FloatImageType::PixelType minPixel(stats->GetMinimum());
      FloatImageType::PixelType maxPixel(stats->GetMaximum());

      if (trimRescaledIntensities > 0.0)
      {
        // REFACTOR: a histogram would be traditional here, but seems
        // over-the-top;
        // I did this because it seemed to me if I knew mean, sigma, max and min,
        // then I know Something about extreme outliers.

        double meanOrig(stats->GetMean());
        double sigmaOrig(stats->GetSigma());

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

        double variationBound((maxPixel - meanOrig) / sigmaOrig);
        double trimBound(variationBound - trimRescaledIntensities);
        if (trimBound > 0.0)
        {
          maxPixel = static_cast<FloatImageType::PixelType>(maxPixel - trimBound * sigmaOrig);
        }
      }

      itk::IntensityWindowingImageFilter<FloatImageType, LOCALImageType>::Pointer remapIntensityFilter =
        itk::IntensityWindowingImageFilter<FloatImageType, LOCALImageType>::New();
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
      itk::CastImageFilter<FloatImageType, LOCALImageType >::Pointer caster = itk::CastImageFilter<FloatImageType, LOCALImageType >::New();
      caster->SetInput(volOrig);
      caster->Update();
      image = caster->GetOutput();
    }
    LOCALImageType::Pointer interpImage;
    if(interpolationMode == "ResampleInPlace")
    {
      using ResampleIPFilterType = itk::ResampleInPlaceImageFilter<SImageType , SImageType>;
      using ResampleIPFilterPointer = ResampleIPFilterType::Pointer;

      using VersorRigid3DTransformType = itk::VersorRigid3DTransform<double>;

      VersorRigid3DTransformType::Pointer result = VersorRigid3DTransformType::New();
      result->SetIdentity();
//      result->SetCenter(orig2msp_img_tfm->GetCenter());
//      result->SetMatrix(orig2msp_img_tfm->GetMatrix());
//      result->SetTranslation(orig2msp_img_tfm->GetTranslation());
      result->Compose(orig2msp_img_tfm);
      ResampleIPFilterPointer resampleIPFilter = ResampleIPFilterType::New();
      resampleIPFilter->SetInputImage(image);
      resampleIPFilter->SetRigidTransform(result);
      resampleIPFilter->Update();
      interpImage = resampleIPFilter->GetOutput();
    }
    else
    {
      // Remember:  the Data is Moving's, the shape is Fixed's.
      interpImage = TransformResample<LOCALImageType, LOCALImageType>(
        image.GetPointer(),
        image.GetPointer(),
        BackgroundFillValue,
        GetInterpolatorFromString<LOCALImageType>(interpolationMode).GetPointer(),
        orig2msp_img_tfm.GetPointer());
    }
    itkUtil::WriteImage<LOCALImageType>(interpImage, resampleMSP);
  }
  if (globalImagedebugLevel > 3)
  {
    const std::string ORIG_ImagePlane(globalResultsDir + "/ORIG_PLANE_" +
                                      itksys::SystemTools::GetFilenameName(inputVolume));
    CreatedebugPlaneImage(internal_image, orig2msp_img_tfm, ORIG_ImagePlane);
  }
  return 0;
}
