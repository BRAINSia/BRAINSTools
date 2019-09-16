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
 * Copyright (c) 2009, Hans J. Johnson
 * All rights reserved.
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
#include "TrimForegroundInDirection.h"
#include "itkLargestForegroundFilledMaskImageFilter.h"
#include <itkImageIteratorWithIndex.h>
#include <itkMinimumMaximumImageCalculator.h>
#include <itkHistogram.h>
#include <itkMath.h>
#include "itkNumberToString.h"

double
FindCenterOfBrainBasedOnTopOfHead(SImageType::Pointer & foreground,
                                  SImageType::Pointer & volOrig,
                                  bool                  maximize,
                                  unsigned int          axis,
                                  double                otsuPercentileThreshold,
                                  unsigned int          closingSize,
                                  double                headSizeLimit,
                                  SImageType::PixelType BackgroundFillValue)
{
  double SI_CenterBasedOnTopOfHead = 0; // This is the return value for the

  // estimated SI location of the center
  // of the brain.

  //   TrimForeground follows these steps:
  //   largest region filled mask for head and neck tissue outline
  //   find maximum Z location of the corner voxels (superior-inferior Z), maxZ
  //   for each voxel in the tissue mask, compute the distance of that voxel
  // PhysicalPoint Z-level to maxZ and store it in a similar code image
  //   make a histogram of the distance map image
  //   count out 2200 CCs to determine a threshold T for too far from the top.
  //   edit the image with the rule, if not in tissue outline or distance map is
  // greater than T, set the source image voxel to zero.

  // ////////////////////////////////////////////////////////////////////////
  //  foreground = FindLargestForgroundFilledMask<SImageType>(volOrig,
  // otsuPercentileThreshold, closingSize);
  using LFFMaskFilterType = itk::LargestForegroundFilledMaskImageFilter<SImageType>;
  LFFMaskFilterType::Pointer LFF = LFFMaskFilterType::New();
  LFF->SetInput(volOrig);
  LFF->SetOtsuPercentileThreshold(otsuPercentileThreshold);
  LFF->SetClosingSize(closingSize);
  LFF->Update();
  foreground = LFF->GetOutput();

  //  foreground will initially hold just a tissue outline (a so-called integer
  // mask),
  //  then it will code for the integer distance to the top of physical space,
  // only within the outline,
  //  then we will convert the zeros and distance codes to input image signal
  // and the assigned background value
  //  in conformity with the assigned volume-distance thresholding from a
  // histogram.

  // ////////////////////////////////////////////////////////////////////////
  //  This will find maximum Superior/Inferior location of the corner voxels.
  double extremum;
  {
    SImageType::SizeType  volOrigSize = volOrig->GetLargestPossibleRegion().GetSize();
    SImageType::IndexType U, V;
    SImageType::PointType limitU, limitV;
    V.Fill(0);
    volOrig->TransformIndexToPhysicalPoint(V, limitV);
    extremum = limitV[axis];
    for (unsigned int i = 0; i < volOrigSize[0]; i += volOrigSize[0] - 1)
    {
      U[0] = i;
      for (unsigned int j = 0; j < volOrigSize[1]; j += volOrigSize[1] - 1)
      {
        U[1] = j;
        for (unsigned int k = 0; k < volOrigSize[2]; k += volOrigSize[2] - 1)
        {
          U[2] = k;
          volOrig->TransformIndexToPhysicalPoint(U, limitU);
          if (maximize ? limitU[axis] > limitV[axis] : limitU[axis] < limitV[axis])
          {
            extremum = limitU[axis];
            limitV = limitU;
            V = U;
          }
        }
      }
    }
    //  Now extremum is the level of the highest subspace plane in the voxel
    // array.
  }
  std::cout << "extremum = " << extremum << std::endl;

  // ////////////////////////////////////////////////////////////////////////
  //  This will produce ForegroundLevel representing where to threshold the head
  // from the neck.
  double ForegroundLevel = 1;
  {
    using SImageIteratorType = itk::ImageRegionIteratorWithIndex<SImageType>;
    SImageIteratorType ItPixel(foreground, foreground->GetLargestPossibleRegion());

    SImageType::PointType PixelPhysicalPoint;
    PixelPhysicalPoint.Fill(0.0);

    ItPixel.Begin();
    for (; !ItPixel.IsAtEnd(); ++ItPixel)
    {
      if (ItPixel.Get() != 0)
      {
        volOrig->TransformIndexToPhysicalPoint(ItPixel.GetIndex(), PixelPhysicalPoint);
        ItPixel.Set(
          static_cast<SImageType::PixelType>(itk::Math::rnd(itk::Math::abs(extremum - PixelPhysicalPoint[axis]))));
      }
      // else, leave the foreground coded zero, not some positive distance from
      // the top.
    }
    // Now foreground holds both the silouette information and the
    // distance-from-extermum information.
  }

  // ////////////////////////////////////////////////////////////////////////
  //  This will populate a histogram to make an increasing volume distribution.
  {
    using HistogramType = itk::Statistics::Histogram<double, 1>;

    using Iterator = itk::ImageRegionIteratorWithIndex<SImageType>;

    double maxval = 0;
    double minval = std::numeric_limits<double>::max();
    using SImageIteratorType = itk::ImageRegionIteratorWithIndex<SImageType>;
    SImageIteratorType imIter(foreground, foreground->GetLargestPossibleRegion());
    while (!imIter.IsAtEnd())
    {
      const double curr_val = imIter.Value();
      if (curr_val > 1) // Need to find min that is greater than zero.
      {
        maxval = (curr_val > maxval) ? curr_val : maxval;
        minval = (curr_val < minval) ? curr_val : minval;
      }
      ++imIter;
    }

    const SImageType::SpacingType origSpacing = volOrig->GetSpacing();
    const double                  voxelSize = origSpacing[0] * origSpacing[1] * origSpacing[2] * 0.001; //
    // Cubic
    // CC.
    std::cout << "voxelSize = " << voxelSize << " cubic CC" << std::endl;
    std::cout << "MIN MAX " << minval << "  " << maxval << std::endl;

    // Make sure that bin width is smaller than the minimum voxel width.  A
    // single layer cannot have more than one row of voxels accumulated in a
    // single bin.
    int numBins = (int)((maxval - minval) / (std::min(origSpacing[0], std::min(origSpacing[1], origSpacing[2]))));

    // Histogram computation
    HistogramType::SizeType size;
    size[0] = numBins;
    std::cout << "Histo values " << minval << " ...  " << maxval << " -> " << numBins << std::endl;


    HistogramType::MeasurementVectorType minValVector, maxValVector;
    minValVector[0] = minval;
    maxValVector[0] = maxval;

    HistogramType::Pointer histogram = HistogramType::New();
    histogram->Initialize(size, minValVector, maxValVector);

    // put each image pixel into the histogram
    HistogramType::MeasurementVectorType measurement;
    Iterator                             iter(foreground, foreground->GetLargestPossibleRegion());
    iter.Begin();
    while (!iter.IsAtEnd())
    {
      const float value = iter.Get();
      measurement[0] = value;
      histogram->IncreaseFrequency(measurement, 1);

      ++iter;
    }

    // ////////////////////////////////////////////////////////////////////////
    //  This will use the histogram to find the desired ForegroundLevel.

    //   SImageType::RegionType imageRegion =
    // foreground->GetLargestPossibleRegion();
    //   int numVoxels = imageRegion.GetSize(0) * imageRegion.GetSize(1) *
    // imageRegion.GetSize(2);
    std::cout << "headSizeLimit = " << headSizeLimit << " CCs" << std::endl;
    double DesiredVolumeToIncludeBeforeClipping = std::numeric_limits<double>::max(); //
    // headSizeLimit
    // is
    // initialized
    // to
    // the
    // smallest
    // possible
    // adult
    // head
    // in
    // cubic
    // cm.

    HistogramType::Iterator           histoIter;
    HistogramType::IndexType          index;
    HistogramType::InstanceIdentifier instance;

    double CummulativeVolume = 0; //  We need to skip over the zero bin.
    bool   exitLoop = false;

    histoIter = histogram->Begin();
    ++histoIter; // Skip the zero bins.

    instance = histoIter.GetInstanceIdentifier();
    index = histogram->GetIndex(instance);
    maxValVector = histogram->GetHistogramMaxFromIndex(index);
    double SupInf_thickness = 0;
    double RLbyAP_area = 0;
    {
      SImageType::PointType physOrigin;
      {
        SImageType::IndexType origin;
        origin[0] = 0;
        origin[1] = 0;
        origin[2] = 0;
        volOrig->TransformIndexToPhysicalPoint(origin, physOrigin);
      }
      SImageType::PointType           physOriginPlusOne;
      itk::ContinuousIndex<double, 3> originPlusOne;
      originPlusOne[0] = volOrig->GetSpacing()[0];
      originPlusOne[1] = volOrig->GetSpacing()[1];
      originPlusOne[2] = volOrig->GetSpacing()[2];
      volOrig->TransformContinuousIndexToPhysicalPoint(originPlusOne, physOriginPlusOne);
      // std::cout << "physOrigin         " << physOrigin        << std::endl;
      // std::cout << "physOriginPlusOne  " << physOriginPlusOne << std::endl;
      const double RL_thickness = itk::Math::abs(physOrigin[0] - physOriginPlusOne[0]) * 0.1;
      const double AP_thickness = itk::Math::abs(physOrigin[1] - physOriginPlusOne[1]) * 0.1;
      SupInf_thickness = itk::Math::abs(physOrigin[2] - physOriginPlusOne[2]) * 0.1; //
      // Convert
      // to
      // cm
      // std::cout << "TEST RL:  " << RL_thickness << " AP " << AP_thickness <<
      // std::endl;
      RLbyAP_area = RL_thickness * AP_thickness; // Convert to cm^2
      if (RLbyAP_area < 1e-5 || SupInf_thickness < 1e-5)
      {
        // std::cout << "  " << SupInf_thickness << std::endl;
        // std::cout << "  " << RL_thickness << std::endl;
        // std::cout << "  " << AP_thickness << std::endl;
        // std::cout << "  " << itk::Math::abs (physOrigin[1]-physOriginPlusOne[1])
        // << std::endl;

        itkGenericExceptionMacro(<< "ERROR:  Can not have zero area, or zero thickness. " << volOrig << std::endl);
      }
    }

    double topOfHeadDistFromExtremeSI = -1;
    double CurrentDistanceFromTopOfHead = 0;
    double largestAreaRadius = 0;
    double MaxCrossSectionalArea = 0.0;
    std::cout << "zero bin count to be skipped = " << histoIter.GetFrequency() << std::endl;
    for (; (histoIter != histogram->End() && !exitLoop); ++histoIter)
    {
      instance = histoIter.GetInstanceIdentifier();
      index = histogram->GetIndex(instance);
      maxValVector = histogram->GetHistogramMaxFromIndex(index);
      minValVector = histogram->GetHistogramMinFromIndex(index);
      if (histoIter.GetFrequency() < 50)
      {
        continue;
      }
      const double CurrentCrossSectionalArea = histoIter.GetFrequency() * RLbyAP_area;
      if (topOfHeadDistFromExtremeSI < 0 && CurrentCrossSectionalArea > 10.0)
      {
        topOfHeadDistFromExtremeSI = maxValVector[0];
      }
      CurrentDistanceFromTopOfHead = (maxValVector[0] - topOfHeadDistFromExtremeSI);
      ForegroundLevel = maxValVector[0];

      if ((CurrentDistanceFromTopOfHead > 70.0) // Require at least 70mm from
                                                // top  of head before
                                                // considering stoping.
          && ((CurrentCrossSectionalArea > MaxCrossSectionalArea) &&
              ((MaxCrossSectionalArea < 10) || (CurrentCrossSectionalArea < MaxCrossSectionalArea * 1.20) //
                                                                                                          // Sometimes
                                                                                                          // histogram
                                                                                                          // bins
                                                                                                          // are
                                                                                                          // filled
                                                                                                          // with
                                                                                                          // 2
                                                                                                          // slices,
                                                                                                          // and
                                                                                                          // that
                                                                                                          // needs
                                                                                                          // to
                                                                                                          // be
                                                                                                          // avoided.
               )))
      {
        MaxCrossSectionalArea = CurrentCrossSectionalArea;
        const double estimated_radius = std::sqrt(MaxCrossSectionalArea / itk::Math::pi); //
        // Estimate
        // the
        // radis
        // of
        // a
        // circle
        // filling
        // this
        // much
        // space
        // Now compute 1.5 times the size of a sphere with this estimated
        // radius.
        constexpr double ScaleFactor = 1.1; // Add 10% for safety
        //  //5+(MaxCrossSectionalArea-200)/100;
        // //Larger brains need more scaling
        const double CurentVolumeBasedOnArea =
          ScaleFactor * (1.33333333333333333 * itk::Math::pi * estimated_radius * estimated_radius * estimated_radius);
        DesiredVolumeToIncludeBeforeClipping = CurentVolumeBasedOnArea;
        // std::cout << "TESTING:  Radius: " << estimated_radius << "
        // DesiredVolume " << DesiredVolumeToIncludeBeforeClipping << std::endl;
      }
      const double CurrentCrossSectionalVolume = histoIter.GetFrequency() * voxelSize;
      CummulativeVolume += CurrentCrossSectionalVolume;
      largestAreaRadius = std::pow(0.75 * itk::Math::one_over_pi * CummulativeVolume, 0.33333333333333333); //
      // Assuming
      // Sphere,
      // what
      // is
      // radius.
      if ((CurrentDistanceFromTopOfHead > 100.0)                         //
                                                                         // Can
                                                                         // not
                                                                         // stop
                                                                         // before
                                                                         // 100
                                                                         // mm
                                                                         // from
                                                                         // top
                                                                         // of
                                                                         // head
                                                                         // are
                                                                         // reached.
          && (CummulativeVolume >= DesiredVolumeToIncludeBeforeClipping) //
                                                                         // Maximum
                                                                         // sustainable
                                                                         // volume
                                                                         // based
                                                                         // on
                                                                         // max
                                                                         // area
                                                                         // of
                                                                         // any
                                                                         // slice.
      )
      {
        std::cout << "VOLUME CRITERIA MET, so exiting. " << CummulativeVolume
                  << " >= " << DesiredVolumeToIncludeBeforeClipping << std::endl;
        exitLoop = true;
      }
      // Now ForegroundLevel represents where to threshold the head from the
      // neck.
    }
    // NOTE:  1 radius was based on some empircal work done by Hans on 100's of
    // data sets.
    SI_CenterBasedOnTopOfHead = extremum - (topOfHeadDistFromExtremeSI + largestAreaRadius * 10.0);
    std::cout << "ForegroundLevel = " << ForegroundLevel << " topOfHeadDistFromExtremeSI " << topOfHeadDistFromExtremeSI
              << " Y_Location_from_Top_of_Head: = " << SI_CenterBasedOnTopOfHead << std::endl;
  }

  // ////////////////////////////////////////////////////////////////////////
  //  This will convert the foreground code image with the rule, foreach voxel,
  //  if not in tissue outline or distance map is greater than T, set the result
  // image voxel to Background;
  //  otherwise set the result image voxel to the source image pixel value.
  {
    using SImageIteratorType = itk::ImageRegionIteratorWithIndex<SImageType>;
    SImageIteratorType ClippedImagePixel(foreground, foreground->GetLargestPossibleRegion());
    SImageIteratorType OriginalImagePixel(volOrig, volOrig->GetLargestPossibleRegion());

    ClippedImagePixel.Begin();
    for (; !ClippedImagePixel.IsAtEnd(); ++ClippedImagePixel)
    {
      if (ClippedImagePixel.Get() != 0)
      {
        if (ClippedImagePixel.Get() <= ForegroundLevel)
        {
          OriginalImagePixel.SetIndex(ClippedImagePixel.GetIndex());
          ClippedImagePixel.Set(OriginalImagePixel.Get());
        }
        else
        {
          ClippedImagePixel.Set(BackgroundFillValue);
        }
      }
      else
      {
        ClippedImagePixel.Set(BackgroundFillValue);
      }
    }
    // Now foreground holds the clipped image.
  }

  return SI_CenterBasedOnTopOfHead;
}
