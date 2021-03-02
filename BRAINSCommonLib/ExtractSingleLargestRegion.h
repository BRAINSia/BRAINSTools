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
#ifndef __ExtractSingleLargestRegion_h
#define __ExtractSingleLargestRegion_h

#include "itkImage.h"

extern itk::Image<unsigned char, 3>::Pointer
ExtractSingleLargestRegionFromMask(const itk::Image<unsigned char, 3>::Pointer & Mask,
                                   const int                                     openingSize,
                                   const int                                     closingSize,
                                   const int                                     safetySize,
                                   const itk::Image<unsigned char, 3>::Pointer & inputLabelImage);

extern itk::Image<unsigned char, 3>::Pointer
ExtractSingleLargestRegion(const unsigned char                           threshold_low,
                           const unsigned char                           threshold_high,
                           const int                                     openingSize,
                           const int                                     closingSize,
                           const int                                     safetySize,
                           const itk::Image<unsigned char, 3>::Pointer & inputLabelImage);

#endif // ExtractSingleLargestRegion_h
