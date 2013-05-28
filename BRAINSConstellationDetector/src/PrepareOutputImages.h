#ifndef __PrepareOutputImages_h__
#define __PrepareOutputImages_h__

#include "GenericTransformImage.h"
#include "landmarksConstellationCommon.h"

namespace itk
{

  extern SImageType::PointType GetOriginalSpaceNamedPoint(
    const LandmarksMapType & landmarks,
    const std::string & NamedPoint);

  extern void PrepareOutputImages(SImageType::Pointer & lOutputResampledImage,
    SImageType::Pointer & lOutputImage,
    SImageType::Pointer & lOutputUntransformedClippedVolume,
    SImageType::ConstPointer lImageToBeResampled,
    VersorTransformType::ConstPointer lVersorTransform,
    VersorTransformType::ConstPointer lInvVersorTransform,
    const LandmarksMapType & inputLmks,
    LandmarksMapType & outputLmks,
    const double lACLowerBound,
    const short int BackgroundFillValue,
    const std::string & lInterpolationMode,
    const bool lCutOutHeadInOutputVolume,
    const double lOtsuPercentileThreshold
  );
}
#endif // __PrepareOutputImages_h__
