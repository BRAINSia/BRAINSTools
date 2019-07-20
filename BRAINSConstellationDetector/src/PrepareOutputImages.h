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
#ifndef __PrepareOutputImages_h__
#define __PrepareOutputImages_h__

#include "GenericTransformImage.h"
#include "landmarksConstellationCommon.h"

namespace itk
{

extern SImageType::PointType
GetNamedPointFromLandmarkList( const LandmarksMapType & landmarks, const std::string & NamedPoint );

extern void
PrepareOutputImages( SImageType::Pointer & lOutputResampledImage, SImageType::Pointer & lOutputImage,
                     SImageType::Pointer &    lOutputUntransformedClippedVolume,
                     SImageType::ConstPointer lImageToBeResampled, VersorTransformType::ConstPointer lVersorTransform,
                     const double lACLowerBound, const short int BackgroundFillValue,
                     const std::string & lInterpolationMode, const bool lCutOutHeadInOutputVolume,
                     const double lOtsuPercentileThreshold );

extern void
ApplyInverseOfTransformToLandmarks( VersorTransformType::ConstPointer lVersorTransform,
                                    const LandmarksMapType & inputLmks, LandmarksMapType & outputLmks );
} // namespace itk
#endif // __PrepareOutputImages_h__
