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
#include <itkImageMomentsCalculator.h>
#include "TrimForegroundInDirection.h"
#include "itkLargestForegroundFilledMaskImageFilter.h"
#include <itkImageIteratorWithIndex.h>
#include <itkMinimumMaximumImageCalculator.h>
#include <itkHistogram.h>
#include <vnl/vnl_math.h>
#include "itkFindCenterOfBrainFilter.h"

// #define USE_DEBUGGIN_IMAGES

void PrintImageInformation(SImageType::Pointer image)
{
  // get information from the image
  SImageType::PointType     origin  = image->GetOrigin();
  SImageType::SpacingType   spacing = image->GetSpacing();
  SImageType::DirectionType direction = image->GetDirection();
  SImageType::IndexType     startIndex = image->GetLargestPossibleRegion().GetIndex();
  SImageType::SizeType      size = image->GetLargestPossibleRegion().GetSize();
  SImageType::PointType     physicalCenter, intensityCenter;
  double                    physicalSize[3] = { 0, 0, 0 };
  SImageType::IndexType     stopIndex;

  SImageType::Pointer infoImage = SImageType::New();    // define an image just
                                                        // to do calculations
  SImageType::RegionType region;
  SImageType::PointType  physicalSpaceMin;
  SImageType::PointType  physicalSpaceMax;
  SImageType::PointType  physicalStartPointSource;
  SImageType::PointType  physicalStopPointSource;

  stopIndex[0] = startIndex[0] + size[0] - 1; // The last valid index
  stopIndex[1] = startIndex[1] + size[1] - 1;
  stopIndex[2] = startIndex[2] + size[2] - 1;
  region.SetSize(size);
  region.SetIndex(startIndex);

  infoImage->SetOrigin(origin); // set the parameters in infoImage
  infoImage->SetSpacing(spacing);
  infoImage->SetDirection(direction);
  infoImage->SetRegions(region);

  // transform startIndex to physical coordinates
  infoImage->TransformIndexToPhysicalPoint(startIndex, physicalStartPointSource);
  infoImage->TransformIndexToPhysicalPoint(stopIndex, physicalStopPointSource);
  physicalSpaceMin = physicalStartPointSource - direction * spacing * 0.5;
  physicalSpaceMax = physicalStopPointSource + direction * spacing * 0.5;

  // physical extent is the difference between max corner and the minimum
  // corner, always positive
  physicalSize[0] = physicalSpaceMax[0] - physicalSpaceMin[0];
  physicalSize[0] = physicalSize[0] > 0 ? physicalSize[0] : 0 - physicalSize[0];
  physicalSize[1] = physicalSpaceMax[1] - physicalSpaceMin[1];
  physicalSize[1] = physicalSize[1] > 0 ? physicalSize[1] : 0 - physicalSize[1];
  physicalSize[2] = physicalSpaceMax[2] - physicalSpaceMin[2];
  physicalSize[2] = physicalSize[2] > 0 ? physicalSize[2] : 0 - physicalSize[2];

  // calculate physical extent as the midpoint between the opposing corners
  physicalCenter[0] = ( physicalSpaceMax[0] + physicalSpaceMin[0] ) / 2;
  physicalCenter[1] = ( physicalSpaceMax[1] + physicalSpaceMin[1] ) / 2;
  physicalCenter[2] = ( physicalSpaceMax[2] + physicalSpaceMin[2] ) / 2;

  // Use the ImageMomentsCalculator to compute the center of mass
  typedef itk::ImageMomentsCalculator<SImageType> Calculator;
  Calculator::Pointer calc = Calculator::New();
  calc->SetImage(image);
  calc->Compute();
  itk::Vector<double, 3> center = calc->GetCenterOfGravity();

  // Output to screen
  std::cout << "IMAGE ORIGIN: " << origin << std::endl;
  std::cout << "IMAGE PHYSICAL EXTENT: [" << physicalSize[0] << ", " << physicalSize[1] << ", " << physicalSize[2]
            << "]" << std::endl;
  std::cout << "IMAGE FIELD OF VIEW: " << physicalSpaceMin << " to " << physicalSpaceMax << std::endl;
  std::cout << "IMAGE CENTER OF PHYSICAL SPACE: " << physicalCenter << std::endl;
  std::cout << "CENTER OF (INTENSITY) MASS: [" << center[0] << ", " << center[1] << ", " << center[2] << "]"
            << std::endl;
}

SImageType::PointType  FindCenterOfBrainBasedOnTopOfHead(SImageType::Pointer & volOrig,
                                                         unsigned int axis,
                                                         double otsuPercentileThreshold,
                                                         unsigned int closingSize,
                                                         double headSizeLowerLimit,
                                                         SImageType::PixelType BackgroundFillValue)
{
  SImageType::Pointer foreground = SImageType::New();

  return TrimForegroundInDirection(foreground,
                                   volOrig,
                                   axis,
                                   otsuPercentileThreshold,
                                   closingSize,
                                   headSizeLowerLimit,
                                   BackgroundFillValue);
}

// ////////////////////////////////////////////////////////////////////////
// ////////////////////////////////////////////////////////////////////////
SImageType::PointType TrimForegroundInDirection(SImageType::Pointer & foreground,  SImageType::Pointer & volOrig,
                                                unsigned int axis,
                                                double otsuPercentileThreshold, unsigned int closingSize,
                                                double headSizeLowerLimit, SImageType::PixelType BackgroundFillValue)
{
  typedef itk::FindCenterOfBrainFilter<SImageType> FindCenterFilter;
  FindCenterFilter::Pointer findCenterFilter = FindCenterFilter::New();
  findCenterFilter->SetInput(volOrig);
  findCenterFilter->SetAxis(axis);
  findCenterFilter->SetOtsuPercentileThreshold(otsuPercentileThreshold);
  findCenterFilter->SetClosingSize(closingSize);
  findCenterFilter->SetHeadSizeLimit(headSizeLowerLimit);
  findCenterFilter->SetBackgroundValue(BackgroundFillValue);

  findCenterFilter->Update();

  foreground = findCenterFilter->GetTrimmedImage();
  SImageType::PointType center = findCenterFilter->GetCenterOfBrain();
  return center;
}
