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
 * Author: Han J. Johnson, Wei Lu, Ali Ghayoor
 * at Psychiatry Imaging Lab,
 * University of Iowa Health Care 2010
 */

#include "itkIO.h"
#include <itkImageMomentsCalculator.h>
#include "itkRecursiveGaussianImageFilter.h"

#include "itkReflectiveCorrelationCenterToImageMetric.h"
#include "landmarksConstellationCommon.h"
#include "TrimForegroundInDirection.h"
#include "itkOrthogonalize3DRotationMatrix.h"

namespace LMC
{
bool debug(false);
bool globalverboseFlag(false);
} // namespace LMC

std::string globalResultsDir("."); // A global variable to define where
                                   // output images are to be placed
int globalImagedebugLevel(1000);   // A global variable to determine the
                                   // level of debugging to perform.

using reflectionFunctorType = Rigid3DCenterReflectorFunctor<itk::PowellOptimizerv4<double>>;

using GaussianFilterType = itk::RecursiveGaussianImageFilter<SImageType, SImageType>;

void
DoMultiQualityReflection(SImageType::Pointer &                  image,
                         RigidTransformType::Pointer &          eyeFixed2msp_lmk_tfm,
                         const int                              qualityLevel,
                         const reflectionFunctorType::Pointer & reflectionFunctor)
{
  // itkUtil::WriteImage<SImageType>(image,"PRE_PYRAMID.nii.gz");
  reflectionFunctor->InitializeImage(image);
  PyramidFilterType::Pointer MyPyramid = MakeThreeLevelPyramid(image.GetPointer());

  SImageType::Pointer EigthImage = MyPyramid->GetOutput(0);
  SImageType::Pointer QuarterImage = MyPyramid->GetOutput(1);
  SImageType::Pointer HalfImage = MyPyramid->GetOutput(2);

  if (qualityLevel >= 0)
  {
    std::cout << "Level 0 Quality Estimates" << std::endl;
    reflectionFunctor->SetDownSampledReferenceImage(EigthImage);
    reflectionFunctor->Initialize();
    reflectionFunctor->Update();
  }
  if (qualityLevel >= 1)
  {
    std::cout << "Level 1 Quality Estimates" << std::endl;
    reflectionFunctor->SetDownSampledReferenceImage(QuarterImage);
    reflectionFunctor->Update();
  }
  if (qualityLevel >= 2)
  {
    std::cout << "Level 2 Quality Estimates" << std::endl;
    reflectionFunctor->SetDownSampledReferenceImage(HalfImage);
    reflectionFunctor->Update();
  }
  if (qualityLevel >= 3)
  {
    std::cout << "Level 3 Quality Estimates" << std::endl;
    reflectionFunctor->SetDownSampledReferenceImage(image);
    reflectionFunctor->Update();
  }
  reflectionFunctor->SetDownSampledReferenceImage(image);
  eyeFixed2msp_lmk_tfm = reflectionFunctor->GetTransformToMSP();
}

void
ComputeMSP(SImageType::Pointer           image,
           RigidTransformType::Pointer & output_transform,
           SImageType::Pointer &         transformedImage,
           const SImageType::PointType & orig_lmk_CenterOfHeadMass,
           const int                     qualityLevel,
           double &                      cc)
{
  if (qualityLevel == -1) // Assume image was pre-aligned outside of the
                          // program
  {
    output_transform = RigidTransformType::New();
    output_transform->SetIdentity();

    itk::ImageDuplicator<SImageType>::Pointer MSP = itk::ImageDuplicator<SImageType>::New();
    MSP->SetInputImage(image);
    MSP->Update();
    transformedImage = MSP->GetOutput();
  }
  else
  {
    reflectionFunctorType::Pointer reflectionFunctor = reflectionFunctorType::New();
    reflectionFunctor->Setorig_lmk_CenterOfHeadMass(orig_lmk_CenterOfHeadMass);

    DoMultiQualityReflection(image, output_transform, qualityLevel, reflectionFunctor);

    transformedImage = reflectionFunctor->GetMSPCenteredImage();
    cc = reflectionFunctor->GetCC();
  }
}

void
ComputeMSP_Easy(SImageType::Pointer           image,
                RigidTransformType::Pointer & eyeFixed2msp_lmk_tfm,
                const SImageType::PointType & orig_lmk_CenterOfHeadMass,
                const int                     qualityLevel)
{
  reflectionFunctorType::Pointer reflectionFunctor = reflectionFunctorType::New();
  reflectionFunctor->Setorig_lmk_CenterOfHeadMass(orig_lmk_CenterOfHeadMass);
  DoMultiQualityReflection(image, eyeFixed2msp_lmk_tfm, qualityLevel, reflectionFunctor);
}

void
CreatedebugPlaneImage(SImageType::Pointer referenceImage, const std::string & debugfilename)
{
  SImageType::PixelType low, high;

  setLowHigh<SImageType>(referenceImage, low, high, 0.01F);

  itk::ImageDuplicator<SImageType>::Pointer MSP = itk::ImageDuplicator<SImageType>::New();
  MSP->SetInputImage(referenceImage);
  MSP->Update();
  SImageType::Pointer           MSPImage = MSP->GetOutput();
  const SImageType::SpacingType imSpacing = MSPImage->GetSpacing();
  SImageType::PointType         CenterOfImage = GetImageCenterPhysicalPoint(MSPImage);

  {
    itk::ImageRegionIteratorWithIndex<SImageType> mspIt(MSPImage, MSPImage->GetLargestPossibleRegion());
    for (; !mspIt.IsAtEnd(); ++mspIt)
    {
      const SImageType::IndexType Index = mspIt.GetIndex();
      SImageType::PointType       Location;
      MSPImage->TransformIndexToPhysicalPoint(Index, Location);
      if (std::abs(Location[0] - CenterOfImage[0]) < imSpacing[0] * 1.00000001)
      {
        mspIt.Set(high);
      }
    }
  }
  itkUtil::WriteImage<SImageType>(MSPImage, debugfilename);
}

SImageType::Pointer
CreatedebugPlaneImage(SImageType::Pointer               referenceImage,
                      const RigidTransformType::Pointer MSPTransform,
                      const std::string &               debugfilename)
{
  SImageType::PixelType low, high;

  setLowHigh<SImageType>(referenceImage, low, high, 0.01F);

  itk::ImageDuplicator<SImageType>::Pointer MSP = itk::ImageDuplicator<SImageType>::New();
  MSP->SetInputImage(referenceImage);
  MSP->Update();
  SImageType::Pointer           MSPImage = MSP->GetOutput();
  const SImageType::SpacingType imSpacing = MSPImage->GetSpacing();
  SImageType::PointType         CenterOfImage = GetImageCenterPhysicalPoint(MSPImage);
  {
    itk::ImageRegionIteratorWithIndex<SImageType> mspIt(MSPImage, MSPImage->GetLargestPossibleRegion());
    for (; !mspIt.IsAtEnd(); ++mspIt)
    {
      const SImageType::IndexType Index = mspIt.GetIndex();
      SImageType::PointType       Location;
      MSPImage->TransformIndexToPhysicalPoint(Index, Location);
      if (std::abs(Location[0] - CenterOfImage[0]) < imSpacing[0] * 1.00000001)
      {
        mspIt.Set(high);
      }
      else
      {
        mspIt.Set(0);
      }
    }
  }
  {
    SImageType::PointType               CrossHairsPoint;
    SImageType::PointType::CoordRepType radius0 = std::abs(3 * imSpacing[0]);
    SImageType::PointType::CoordRepType radius1 = std::abs(3 * imSpacing[1]);
    SImageType::PointType::CoordRepType radius2 = std::abs(3 * imSpacing[2]);
    for (SImageType::PointType::CoordRepType k = CenterOfImage[2] - radius2; k < CenterOfImage[2] + radius2;
         k += imSpacing[2])
    {
      for (SImageType::PointType::CoordRepType j = CenterOfImage[1] - radius1; j < CenterOfImage[1] + radius1;
           j += imSpacing[1])
      {
        for (SImageType::PointType::CoordRepType i = CenterOfImage[0] - radius0; i < CenterOfImage[0] + radius0;
             i += imSpacing[0])
        {
          CrossHairsPoint[0] = i;
          CrossHairsPoint[1] = j;
          CrossHairsPoint[2] = k;
          SImageType::IndexType closestIndex;
          const bool            isInside = MSPImage->TransformPhysicalPointToIndex(CrossHairsPoint, closestIndex);
          if (isInside)
          {
            MSPImage->SetPixel(closestIndex, high);
          }
        }
      }
    }
  }

  RigidTransformType::Pointer invMSPTransform = RigidTransformType::New();
  MSPTransform->GetInverse(invMSPTransform);
  SImageType::Pointer RotatedPlane =
    TransformResample<SImageType, SImageType>(MSPImage.GetPointer(),
                                              MSPImage.GetPointer(),
                                              /* DEFAULT TO ZERO */ 0,
                                              GetInterpolatorFromString<SImageType>("Linear").GetPointer(),
                                              invMSPTransform.GetPointer());

  itk::ImageDuplicator<SImageType>::Pointer duplicator = itk::ImageDuplicator<SImageType>::New();
  duplicator->SetInputImage(referenceImage);
  duplicator->Update();
  SImageType::Pointer RasterImage = duplicator->GetOutput();

  itk::ImageRegionIteratorWithIndex<SImageType> rasterIt(RasterImage, RasterImage->GetLargestPossibleRegion());
  itk::ImageRegionIteratorWithIndex<SImageType> rplaneIt(RotatedPlane, RotatedPlane->GetLargestPossibleRegion());
  for (rplaneIt.GoToBegin(); !rplaneIt.IsAtEnd() && !rasterIt.IsAtEnd(); ++rplaneIt, ++rasterIt)
  {
    if (rasterIt.Get() > high)
    {
      rasterIt.Set(high);
    }
    if (rplaneIt.Get() > high * 0.5)
    {
      SImageType::PixelType p = static_cast<SImageType::PixelType>((2.0 * rplaneIt.Get() + rasterIt.Get()) * 0.3333333);
      rasterIt.Set(p);
    }
  }
  itkUtil::WriteImage<SImageType>(RasterImage, debugfilename);
  return RasterImage;
}


PyramidFilterType::Pointer
MakeThreeLevelPyramid(SImageType::Pointer refImage)
{
  PyramidFilterType::ScheduleType pyramidSchedule;

  PyramidFilterType::Pointer MyPyramid = PyramidFilterType::New();

  MyPyramid->SetInput(refImage);
  MyPyramid->SetNumberOfLevels(3);
  pyramidSchedule.SetSize(3, 3);
  // Attempt to set a schedule so that the top of the pyramid
  // has images of about 8mm, and the next level has resolutions about 4mm
  // isotropic voxels
  // these are sizes found to work well for estimating MSP without making the
  // image too small.
  SImageType::SpacingType refImageSpacing = refImage->GetSpacing();
  for (unsigned int c = 0; c < pyramidSchedule.cols(); ++c)
  {
    // about 8mm
    pyramidSchedule[0][c] = static_cast<unsigned int>(2 * round(4.0 / refImageSpacing[c]));
    // about 4mm
    pyramidSchedule[1][c] = static_cast<unsigned int>(2 * round(2.0 / refImageSpacing[c]));
    // about 2mm
    pyramidSchedule[2][c] = static_cast<unsigned int>(2 * round(1.0 / refImageSpacing[c]));
  }
  MyPyramid->SetSchedule(pyramidSchedule);
  MyPyramid->Update();
  return MyPyramid;
}

PyramidFilterType::Pointer
MakeOneLevelPyramid(SImageType::Pointer refImage)
{
  PyramidFilterType::ScheduleType pyramidSchedule;

  PyramidFilterType::Pointer MyPyramid = PyramidFilterType::New();

  MyPyramid->SetInput(refImage);
  MyPyramid->SetNumberOfLevels(1);
  pyramidSchedule.SetSize(1, 3);

  SImageType::SpacingType refImageSpacing = refImage->GetSpacing();
  for (unsigned int c = 0; c < pyramidSchedule.cols(); ++c)
  {
    // about 8mm
    pyramidSchedule[0][c] = static_cast<unsigned int>(2 * round(4.0 / refImageSpacing[c]));
  }
  MyPyramid->SetSchedule(pyramidSchedule);
  MyPyramid->Update();
  return MyPyramid;
}

// ////////////////////////////////////////////////////////////////////////
//  This is a lightweight wrapper for FindCenterOfBrainBasedOnTopOfHead, which
// evolved from TrimForegroundInDirection.
SImageType::PointType
GetCenterOfHeadMass(SImageType::Pointer volume)
{
  // NOTE:   This is a bit of a hack to deal with vastly different sized heads,
  // note:   or where the field of view is not centered in the brain
  // note:   or where there is a massive amount of neck.

  // Choose the Inferior/Superior based on maximum dirction cosign

  // Find center of the image space which is where the MSP has been placed.
  constexpr unsigned int ISdirectionIndex = 2;
  SImageType::PointType  CenterOfMass = FindCenterOfBrainBasedOnTopOfHead(volume, ISdirectionIndex, 0.01, 7, 700, 0);

  // IN AN LPS SYSTEM, the INFERIOR/SUPERIOR should be the center of physical
  // space, rather than the
  // center of mass because the INFERIOR/SUPERIOR direction depends so much on
  // the size of the neck.

  return CenterOfMass;
}

//
//
//
// ////////////////////////////////////////////////////////////////////////////////////////////////
RigidTransformType::Pointer
computeTmspFromPoints(SImageType::PointType RP,
                      SImageType::PointType AC,
                      SImageType::PointType PC,
                      SImageType::PointType DesiredCenter)
{
  // a variable to store correlation coefficient values
  SImageType::PointType::VectorType ACPC = PC - AC;

  ACPC.Normalize();
  SImageType::PointType::VectorType ACRP = RP - AC;

  vnl_vector_fixed<double, 3> NormalToMSPPlane = vnl_cross_3d(ACPC.GetVnlVector(), ACRP.GetVnlVector());
  // --std::cout << ACPC << "=" << PC << " - " << AC << std::endl;
  // --std::cout << ACRP << "=" << RP << " - " << AC << std::endl;
  // --std::cout << NormalToMSPPlane << std::endl;
  SImageType::PointType::VectorType NormalToMSPPlaneVector;
  NormalToMSPPlaneVector[0] = NormalToMSPPlane[0];
  NormalToMSPPlaneVector[1] = NormalToMSPPlane[1];
  NormalToMSPPlaneVector[2] = NormalToMSPPlane[2];
  NormalToMSPPlaneVector.Normalize();
  if (NormalToMSPPlaneVector[0] < 0)
  {
    NormalToMSPPlaneVector *= -1.0;
  }

  // --std::cout << "Norm: " << NormalToMSPPlaneVector << "     NORMACPC: " <<
  // ACPC << std::endl;

  // double
  // PlaneNormalBank=-std::atan2(NormalToMSPPlaneVector[2],NormalToMSPPlaneVector[0]);
  //  //Rotate the "Y" (i.e. Anterior to Posterior Axis)
  // double PlaneNormalHeading=-std::acos(NormalToMSPPlaneVector[0]);
  //               //Rotate the "Z" (i.e. Inferior to Superior Axis)

  double PlaneNormalBank = -std::atan2(NormalToMSPPlaneVector[2], NormalToMSPPlaneVector[0]);
  // Rotate the "Y" (i.e. Anterior to Posterior Axis)
  double PlaneNormalHeading = std::sin(NormalToMSPPlaneVector[1]);
  // Rotate the "Z" (i.e. Inferior to Superior Axis)
  double PlaneNormalAttitude = std::sin(ACPC[2]);

  SImageType::PointType::VectorType CenterOffset = AC.GetVectorFromOrigin() - DesiredCenter.GetVectorFromOrigin();

  RigidTransformType::Pointer AlignMSPTransform = RigidTransformType::New();
  AlignMSPTransform->SetCenter(DesiredCenter);
  AlignMSPTransform->SetRotation(PlaneNormalAttitude, PlaneNormalBank, PlaneNormalHeading);
  {
    // Clean up the rotation to make it orthogonal:
    const itk::Matrix<double, 3, 3> & CleanedOrthogonalized =
      itk::Orthogonalize3DRotationMatrix(AlignMSPTransform->GetMatrix());
    AlignMSPTransform->SetMatrix(CleanedOrthogonalized);
  }
  AlignMSPTransform->SetTranslation(CenterOffset);
  if (LMC::globalverboseFlag)
  {
    std::cout << "\n============================" << CenterOffset << " \n" << AlignMSPTransform << std::endl;
    std::cout << " \n" << AlignMSPTransform->GetTranslation() << std::endl;
    std::cout << "============================\n" << std::endl;
  }
  return AlignMSPTransform;
}


//
//
// ////////////////////////////////////////////////////////////////////////////////////////////////

// Use the law of cosines to determine the angle between the RPtoPC vector and
// the RPtoAC vector.
void
decomposeRPAC(const SImageType::PointType & RP,
              const SImageType::PointType & PC,
              const SImageType::PointType & AC,
              double * const                RPPC_to_RPAC_angle,
              double * const                RPAC_over_RPPC)
{
  const double RPtoPC = std::sqrt(((RP[1] - PC[1]) * (RP[1] - PC[1]) + (RP[2] - PC[2]) * (RP[2] - PC[2])));
  const double RPtoAC = std::sqrt(((RP[1] - AC[1]) * (RP[1] - AC[1]) + (RP[2] - AC[2]) * (RP[2] - AC[2])));
  const double PCtoAC = std::sqrt(((PC[1] - AC[1]) * (PC[1] - AC[1]) + (PC[2] - AC[2]) * (PC[2] - AC[2])));

  const double cos_gamma = (RPtoPC * RPtoPC + RPtoAC * RPtoAC - PCtoAC * PCtoAC) / (2.0 * RPtoPC * RPtoAC);
  const double gamma = std::acos(cos_gamma);

  // const double sin_gamma = std::sqrt( 1.0 - cos_gamma*cos_gamma );

  *RPPC_to_RPAC_angle = gamma;       // The angle between the RP-PC vector and
                                     // the RP-AC vector
  *RPAC_over_RPPC = RPtoAC / RPtoPC; // The ratio of lengths between the RP-PC
                                     // vector and the RP-AC vector
                                     // std::cout << "----------TEST ANGLE " <<
                                     //  *RPPC_to_RPAC_angle *180.0/itk::Math::pi
                                     // << "  RATIO " << *RPAC_over_RPPC <<
                                     // std::endl;
}


void
extractArray(LinearInterpolatorType::Pointer                                imInterp,
             const SImageType::PointType &                                  CenterPoint,
             const landmarksConstellationModelIO::IndexLocationVectorType & model,
             std::vector<float> &                                           result_array)
{
  int q = 0;
  for (landmarksConstellationModelIO::IndexLocationVectorType::const_iterator it = model.begin(); it != model.end();
       ++it, ++q)
  {
    const SImageType::PointType & point = CenterPoint + *it;
    if (imInterp->IsInsideBuffer(point))
    {
      result_array[q] = static_cast<float>(imInterp->Evaluate(point));
    }
    else
    {
      result_array[q] = 0.0F;
    }
  }
}

unsigned char
ShortToUChar(short in, short min, short max)
{
  double x(in);

  x -= static_cast<double>(min);
  double divisor = max - min;
  x /= divisor;
  x *= 255.0;
  return static_cast<unsigned char>(x);
}

SImageType::Pointer
CreateTestCenteredRotatedImage2(const RigidTransformType::Pointer ACPC_MSP_AlignedTransform,
                                /* const
                                  SImageType::PointType
                                  finalPoint_MSP, */
                                const SImageType::PointType         itkNotUsed(PreMSP_Point),
                                /*const*/ SImageType::Pointer &     image,
                                const RigidTransformType::Pointer & Point_Rotate)
{
  // ////// Compose test rotation with the translated ACPC alignment
  RigidTransformType::Pointer Point_Centered_TestRotated = RigidTransformType::New();

  Point_Centered_TestRotated->SetFixedParameters(ACPC_MSP_AlignedTransform->GetFixedParameters());
  Point_Centered_TestRotated->SetParameters(ACPC_MSP_AlignedTransform->GetParameters());
  // INFO:  Perhaps change sign of Point_Translate, and the remove the true flag.
  Point_Centered_TestRotated->Compose(Point_Rotate, true);


  SImageType::Pointer image_Point_TestRotated =
    TransformResample<SImageType, SImageType>(image.GetPointer(),
                                              image.GetPointer(),
                                              /* DEFAULT TO ZERO */ 0,
                                              GetInterpolatorFromString<SImageType>("Linear").GetPointer(),
                                              Point_Centered_TestRotated.GetPointer());

  RigidTransformType::Pointer invPoint_Centered_TestRotated = RigidTransformType::New();
  Point_Centered_TestRotated->GetInverse(invPoint_Centered_TestRotated);
  return image_Point_TestRotated;
}

// using PointType = landmarksDataSet::PointType;

void
MakeLabelImage(SImageType::Pointer           in,
               const SImageType::PointType & RP,
               const SImageType::PointType & AC,
               const SImageType::PointType & PC,
               const SImageType::PointType & VN4,
               const std::string &           fname)
{
  SImageType::Pointer maskImage = SImageType::New();

  maskImage->CopyInformation(in);
  maskImage->SetRegions(in->GetLargestPossibleRegion());
  maskImage->Allocate();
  maskImage->FillBuffer(0);
  {
    SImageType::IndexType Index;
    maskImage->TransformPhysicalPointToIndex(RP, Index);
    maskImage->SetPixel(Index, 250);
    maskImage->TransformPhysicalPointToIndex(AC, Index);
    maskImage->SetPixel(Index, 255);
    maskImage->TransformPhysicalPointToIndex(PC, Index);
    maskImage->SetPixel(Index, 245);
    maskImage->TransformPhysicalPointToIndex(VN4, Index);
    maskImage->SetPixel(Index, 240);
  }
  itkUtil::WriteImage<SImageType>(maskImage, fname);
}

SImageType::PointType
GetImageCenterPhysicalPoint(SImageType::Pointer & image)
{
  const SImageType::SizeType imageOverallSize = image->GetLargestPossibleRegion().GetSize();

  itk::ContinuousIndex<double, 3> centerIndex;

  for (size_t q = 0; q < SImageType::ImageDimension; ++q)
  {
    centerIndex[q] = 0.5 * (imageOverallSize[q] - 1);
  }
  SImageType::PointType centerLocation;
  image->TransformContinuousIndexToPhysicalPoint(centerIndex, centerLocation);
  return centerLocation;
}

// INFO:  Need to make size and spacing command line arguments that are sent into
// the ResampleToIsotropicImage filter.
// INFO:  Currently hard-coded to 1.0^3 isotropic voxels in a 256^3 matrix
// --outputVolumeSize 256,256,256 --outputVolumeSpacing 1.0,1.0,1.0
SImageType::Pointer
MakeIsoTropicReferenceImage()
{
  SImageType::DirectionType Ident;

  Ident.SetIdentity();
  SImageType::PointType Origin;
  Origin.Fill(-127.5);
  SImageType::SpacingType Spacing;
  Spacing.Fill(1.0);

  SImageType::SizeType Size;
  Size.Fill(256);
  SImageType::RegionType Region;
  Region.SetSize(Size);

  SImageType::Pointer isotropicReferenceVolume = SImageType::New();
  isotropicReferenceVolume->SetDirection(Ident);
  isotropicReferenceVolume->SetOrigin(Origin);
  isotropicReferenceVolume->SetSpacing(Spacing);
  isotropicReferenceVolume->SetRegions(Region);
  // NOTE:  No need to allocate the image memory, it is just used as a template.
  return isotropicReferenceVolume;
}