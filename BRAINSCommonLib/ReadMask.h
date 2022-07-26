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
#ifndef ReadMask_h_
#define ReadMask_h_

#include "itkIO.h"
#include "itkImageMaskSpatialObject.h"
#include "itkBinaryThresholdImageFilter.h"

template <typename MaskType, unsigned VDimension>
typename MaskType::Pointer
ReadImageMask(const std::string & filename, typename itk::ImageBase<VDimension> * /*referenceImage*/)
{
  using MaskPixelType = unsigned char;
  using SpatialObjectMaskPixelType = itk::Image<MaskPixelType, 3 >;

  // convert mask image to mask
  using ReadImageMaskSpatialObjectType = itk::ImageMaskSpatialObject<SpatialObjectMaskPixelType::ImageDimension>;
  typename ReadImageMaskSpatialObjectType::Pointer mask = ReadImageMaskSpatialObjectType::New();
  {
    // Allow reading double precision images, and rely on standard typecasting to map zero to zero, and non-zero to non-zero.
    using ReadMaskImageType = itk::Image<double, 3>;
    typename ReadMaskImageType::Pointer OrientedMaskImage = itkUtil::ReadImage<ReadMaskImageType>(filename);
    using BinaryThresholdFilterType = itk::BinaryThresholdImageFilter<ReadMaskImageType, SpatialObjectMaskPixelType>;
    BinaryThresholdFilterType::Pointer btf = BinaryThresholdFilterType::New();
    btf->SetInput(OrientedMaskImage);
    //NOTE: Identifying the background as Zero, so that all non-zero is foreground.
    btf->SetLowerThreshold(0);
    btf->SetUpperThreshold(0);
    btf->SetInsideValue(0);
    btf->SetOutsideValue(1);
    btf->Update();
    mask->SetImage(btf->GetOutput());
  }

  mask->Update(); // Replaced old ComputeObjectToWorldTransform with new Update()
  // return pointer to mask
  typename MaskType::Pointer p = dynamic_cast<MaskType *>(mask.GetPointer());
  if (p.IsNull())
  {
    itkGenericExceptionMacro(<< "Failed conversion to Mask");
  }
  return p;
}

#endif // LoadMask_h
