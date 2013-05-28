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
    const double lACLowerBound,
    const short int BackgroundFillValue,
    const std::string & lInterpolationMode,
    const bool lCutOutHeadInOutputVolume,
    const double lOtsuPercentileThreshold
  );

  extern void PrepareOutputLandmarks(
    VersorTransformType::ConstPointer lVersorTransform,
    const LandmarksMapType & inputLmks,
    LandmarksMapType & outputLmks
    );
}
#endif // __PrepareOutputImages_h__
