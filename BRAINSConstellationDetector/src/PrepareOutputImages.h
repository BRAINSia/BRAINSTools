#ifndef __PrepareOutputImages_h__

#include "BRAINSConstellationDetector2.h"
#include "GenericTransformImage.h"
#include "itkOrthogonalize3DRotationMatrix.h"


namespace itk
{

extern void PrepareOutputImages(SImageType::Pointer & lOutputResampledImage,
    SImageType::Pointer & lOutputImage,
    SImageType::Pointer & lOutputUntransformedClippedVolume,
    SImageType::Pointer & lImageToBeResampled,
    VersorTransformType::Pointer & lInvVersorTransform,
    const double lACLowerBound,
    const double BackgroundFillValue,
    const double PhysicalLowerBound,
    const bool lCutOutHeadInOutputVolume,
    const double lOtsuPercentileThreshold
  );
}
#endif // __PrepareOutputImages_h__
