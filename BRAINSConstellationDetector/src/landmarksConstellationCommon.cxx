/*
 * Author: Han J. Johnson, Wei Lu
 * at Psychiatry Imaging Lab,
 * University of Iowa Health Care 2010
 */

#include "itkIO.h"
#include <itkImageMomentsCalculator.h>

#include "itkReflectiveCorrelationCenterToImageMetric.h"
#include "landmarksConstellationCommon.h"
#include "TrimForegroundInDirection.h"
#include "GenericTransformImage.h"
#include "itkOrthogonalize3DRotationMatrix.h"

typedef itk::ImageMomentsCalculator<SImageType> momentsCalculatorType;

// TODO:  This should really be a singleton.
RandomGeneratorPointer _RandomGenerator;
const unsigned int     MAX_ROTATIONS_TESTED = 90;
const unsigned int     MAXITER = 5000;
const unsigned int     DEL = 3;
const unsigned int     YES = 1;
const unsigned int     NO = 0;
const unsigned int     SMAX = 50;
namespace LMC
{
bool debug(false);
bool globalverboseFlag(false);
}

std::string globalResultsDir(".");       // A global variable to define where
                                         // output images are to be placed
int globalImagedebugLevel(1000);         // A global variable to determine the
                                         // level of debugging to perform.

void ComputeMSP(SImageType::Pointer image,
                RigidTransformType::Pointer & Tmsp,
                SImageType::Pointer & transformedImage,
                const SImageType::PointType & centerOfHeadMass,
                const int qualityLevel,
                double & cc)
{
  if( qualityLevel == -1 )  // Assume image was pre-aligned outside of the
                            // program
    {
    Tmsp = RigidTransformType::New();
    Tmsp->SetIdentity();

    itk::ImageDuplicator<SImageType>::Pointer MSP = itk::ImageDuplicator<SImageType>::New();
    MSP->SetInputImage(image);
    MSP->Update();
    transformedImage = MSP->GetModifiableOutput();
    }
  else
    {
    // itkUtil::WriteImage<SImageType>(image,"PRE_PYRAMID.nii.gz");
    PyramidFilterType::Pointer    MyPyramid = MakeThreeLevelPyramid(image.GetPointer() );
    SImageType::Pointer           EigthImage = MyPyramid->GetOutput(0);
    SImageType::Pointer           QuarterImage = MyPyramid->GetOutput(1);
    SImageType::Pointer           HalfImage = MyPyramid->GetOutput(2);
    Rigid3DCenterReflectorFunctor reflectionFunctor;

    reflectionFunctor.SetCenterOfHeadMass(centerOfHeadMass);
    reflectionFunctor.Initialize(image);
    if( qualityLevel >= 0 )
      {
      std::cout << "Level 0 Quality Estimates" << std::endl;
      reflectionFunctor.SetDownSampledReferenceImage(EigthImage);
      reflectionFunctor.QuickSampleParameterSpace();
      reflectionFunctor.Update();
      }
    if( qualityLevel >= 1 )
      {
      std::cout << "Level 1 Quality Estimates" << std::endl;
      reflectionFunctor.SetDownSampledReferenceImage(QuarterImage);
      reflectionFunctor.Update();
      }
    if( qualityLevel >= 2 )
      {
      std::cout << "Level 2 Quality Estimates" << std::endl;
      reflectionFunctor.SetDownSampledReferenceImage(HalfImage);
      reflectionFunctor.Update();
      }
    if( qualityLevel >= 3 )
      {
      std::cout << "Level 3 Quality Estimates" << std::endl;
      reflectionFunctor.SetDownSampledReferenceImage(image);
      reflectionFunctor.Update();
      }
    reflectionFunctor.SetDownSampledReferenceImage(image);
    Tmsp = reflectionFunctor.GetTransformToMSP();
    transformedImage = reflectionFunctor.GetMSPCenteredImage();
    cc = reflectionFunctor.GetCC();
    }
}

void ComputeMSP(SImageType::Pointer image, RigidTransformType::Pointer & Tmsp, const int qualityLevel)
{
  // itkUtil::WriteImage<SImageType>(image,"PRE_PYRAMID.nii.gz");
  PyramidFilterType::Pointer    MyPyramid = MakeThreeLevelPyramid(image);
  SImageType::Pointer           EigthImage = MyPyramid->GetOutput(0);
  SImageType::Pointer           QuarterImage = MyPyramid->GetOutput(1);
  SImageType::Pointer           HalfImage = MyPyramid->GetOutput(2);
  Rigid3DCenterReflectorFunctor reflectionFunctor;

  reflectionFunctor.Initialize(image);
  if( qualityLevel >= 0 )
    {
    std::cout << "Level 0 Quality Estimates" << std::endl;
    reflectionFunctor.SetDownSampledReferenceImage(EigthImage);
    reflectionFunctor.QuickSampleParameterSpace();
    reflectionFunctor.Update();
    }
  if( qualityLevel >= 1 )
    {
    std::cout << "Level 1 Quality Estimates" << std::endl;
    reflectionFunctor.SetDownSampledReferenceImage(QuarterImage);
    reflectionFunctor.Update();
    }
  if( qualityLevel >= 2 )
    {
    std::cout << "Level 2 Quality Estimates" << std::endl;
    reflectionFunctor.SetDownSampledReferenceImage(HalfImage);
    reflectionFunctor.Update();
    }
  if( qualityLevel >= 3 )
    {
    std::cout << "Level 3 Quality Estimates" << std::endl;
    reflectionFunctor.SetDownSampledReferenceImage(image);
    reflectionFunctor.Update();
    }
  reflectionFunctor.SetDownSampledReferenceImage(image);
  Tmsp = reflectionFunctor.GetTransformToMSP();
}

void CreatedebugPlaneImage(SImageType::Pointer referenceImage, const std::string & debugfilename)
{
  SImageType::PixelType low, high;

  setLowHigh<SImageType>(referenceImage, low, high, 0.01F);

  itk::ImageDuplicator<SImageType>::Pointer MSP = itk::ImageDuplicator<SImageType>::New();
  MSP->SetInputImage(referenceImage);
  MSP->Update();
  SImageType::Pointer           MSPImage = MSP->GetModifiableOutput();
  const SImageType::SpacingType imSpacing = MSPImage->GetSpacing();
  SImageType::PointType         CenterOfImage = GetImageCenterPhysicalPoint(MSPImage);

    {
    itk::ImageRegionIteratorWithIndex<SImageType> mspIt( MSPImage, MSPImage->GetLargestPossibleRegion() );
    for( ; !mspIt.IsAtEnd(); ++mspIt )
      {
      const SImageType::IndexType Index = mspIt.GetIndex();
      SImageType::PointType       Location;
      MSPImage->TransformIndexToPhysicalPoint(Index, Location);
      if( vcl_abs(Location[0] - CenterOfImage[0]) < imSpacing[0] * 1.00000001 )
        {
        mspIt.Set(high);
        }
      }
    }
  itkUtil::WriteImage<SImageType>(MSPImage, debugfilename);
}

// TODO:  Should be ::ConstPointer
SImageType::Pointer CreatedebugPlaneImage(SImageType::Pointer referenceImage,
                                          const RigidTransformType::Pointer MSPTransform,
                                          const std::string & debugfilename)
{
  SImageType::PixelType low, high;

  setLowHigh<SImageType>(referenceImage, low, high, 0.01F);

  itk::ImageDuplicator<SImageType>::Pointer MSP = itk::ImageDuplicator<SImageType>::New();
  MSP->SetInputImage(referenceImage);
  MSP->Update();
  SImageType::Pointer           MSPImage = MSP->GetModifiableOutput();
  const SImageType::SpacingType imSpacing = MSPImage->GetSpacing();
  SImageType::PointType CenterOfImage = GetImageCenterPhysicalPoint(MSPImage);
    {
    itk::ImageRegionIteratorWithIndex<SImageType> mspIt( MSPImage, MSPImage->GetLargestPossibleRegion() );
    for( ; !mspIt.IsAtEnd(); ++mspIt )
      {
      const SImageType::IndexType Index = mspIt.GetIndex();
      SImageType::PointType       Location;
      MSPImage->TransformIndexToPhysicalPoint(Index, Location);
      if( vcl_abs(Location[0] - CenterOfImage[0]) < imSpacing[0] * 1.00000001 )
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
    SImageType::PointType::CoordRepType radius0 = vcl_abs(3 * imSpacing[0]);
    SImageType::PointType::CoordRepType radius1 = vcl_abs(3 * imSpacing[1]);
    SImageType::PointType::CoordRepType radius2 = vcl_abs(3 * imSpacing[2]);
    for( SImageType::PointType::CoordRepType k = CenterOfImage[2] - radius2;
         k < CenterOfImage[2] + radius2;
         k += imSpacing[2] )
      {
      for( SImageType::PointType::CoordRepType j = CenterOfImage[1] - radius1;
           j < CenterOfImage[1] + radius1;
           j += imSpacing[1] )
        {
        for( SImageType::PointType::CoordRepType i = CenterOfImage[0] - radius0;
             i < CenterOfImage[0] + radius0;
             i += imSpacing[0] )
          {
          CrossHairsPoint[0] = i; CrossHairsPoint[1] = j; CrossHairsPoint[2] = k;
          SImageType::IndexType closestIndex;
          const bool            isInside = MSPImage->TransformPhysicalPointToIndex(CrossHairsPoint, closestIndex);
          if( isInside )
            {
            MSPImage->SetPixel(closestIndex, high);
            }
          }
        }
      }
    }

  RigidTransformType::Pointer invMSPTransform = RigidTransformType::New();
  MSPTransform->GetInverse(invMSPTransform);
#if 1
  SImageType::Pointer RotatedPlane = TransformResample<SImageType, SImageType>(
      MSPImage, MSPImage, /* DEFAULT TO ZERO */ 0,
      GetInterpolatorFromString<SImageType>("Linear"), invMSPTransform.GetPointer() );
#else
  SImageType::Pointer RotatedPlane = GenericImageResample<SImageType, SImageType, LinearInterpolatorType>(
      MSPImage,
      invMSPTransform,
      /* DEFAULT TO ZERO */
      0);
#endif

  itk::ImageDuplicator<SImageType>::Pointer duplicator = itk::ImageDuplicator<SImageType>::New();
  duplicator->SetInputImage(referenceImage);
  duplicator->Update();
  SImageType::Pointer RasterImage = duplicator->GetModifiableOutput();

  itk::ImageRegionIteratorWithIndex<SImageType> rasterIt( RasterImage, RasterImage->GetLargestPossibleRegion() );
  itk::ImageRegionIteratorWithIndex<SImageType> rplaneIt( RotatedPlane, RotatedPlane->GetLargestPossibleRegion() );
  for( rplaneIt.GoToBegin(); !rplaneIt.IsAtEnd() && !rasterIt.IsAtEnd(); ++rplaneIt, ++rasterIt )
    {
    if( rasterIt.Get() > high )
      {
      rasterIt.Set(high);
      }
    if( rplaneIt.Get() > high * 0.5 )
      {
      SImageType::PixelType p =
        static_cast<SImageType::PixelType>
        ( ( 2.0 * rplaneIt.Get() + rasterIt.Get() ) * 0.3333333 );
      rasterIt.Set(p);
      }
    }
  itkUtil::WriteImage<SImageType>(RasterImage, debugfilename);
  return RasterImage;
}

vnl_vector<double> convertToReadable(const vnl_vector<double> & input)
{
  vnl_vector<double> temp;
  temp.set_size(3);
  temp[0] = input[0] * 180.0 / vnl_math::pi;
  temp[1] = input[1] * 180.0 / vnl_math::pi;
  temp[2] = input[2];
  return temp;
}

/** this conversion uses conventions as described on page:
 *   http://www.euclideanspace.com/maths/geometry/rotations/euler/index.htm
 *   Coordinate System: right hand
 *   Positive angle: right hand
 *   Order of euler angles: initialHeadingAngle first, then initialAttitudeAngle, then initialBankAngle
 *   matrix row column ordering:
 *   [m00 m01 m02]
 *   [m10 m11 m12]
 *   [m20 m21 m22]*/
// With respect to an LPS system of the head,
// initialHeadingAngle is where you are looking left to right, (Rotation of Z
// axis)
// attiude is where you are looking from floor to ceiling (Rotation of X axis)
// initialBankAngle is how tilted your head is  (Rotation of Y axis)
void ComputeEulerAnglesFromRotationMatrix(const itk::Matrix<double, 3, 3> &  m,
                                          double & initialAttitudeAngle,
                                          double & initialBankAngle,
                                          double & initialHeadingAngle)
{
  // Assuming the angles are in radians.
  if( m[1][0] > 0.998 )  // singularity at north pole
    {
    initialAttitudeAngle = 0;
    initialBankAngle = vcl_atan2(m[0][2], m[2][2]);
    initialHeadingAngle = vnl_math::pi_over_2;
    return;
    }
  if( m[1][0] < -0.998 )  // singularity at south pole
    {
    initialAttitudeAngle = 0;
    initialBankAngle = vcl_atan2(m[0][2], m[2][2]);
    initialHeadingAngle = -vnl_math::pi_over_2;
    return;
    }
  initialAttitudeAngle = vcl_atan2(-m[1][2], m[1][1]);
  initialBankAngle = vcl_atan2(-m[2][0], m[0][0]);
  initialHeadingAngle = vcl_asin(m[1][0]);
}


itk::Versor<double> CreateRotationVersorFromAngles(const double alpha, const double beta, const double gamma)
{
  // http://en.wikipedia.org/wiki/Conversion_between_quaternions_and_Euler_angles
  // psi = alpha is rotate the X axis -- Attitude
  // theta= beta is rotate the Y axis  -- Bank
  // phi=  gamma is rotate the Z axis -- Heading
  const double cha = vcl_cos(alpha * 0.5);
  const double chb = vcl_cos(beta * 0.5);
  const double chg = vcl_cos(gamma * 0.5);
  const double sha = vcl_sin(alpha * 0.5);
  const double shb = vcl_sin(beta * 0.5);
  const double shg = vcl_sin(gamma * 0.5);

  vnl_vector_fixed<double, 4> q;
  q[0] = cha * chb * chg + sha * shb * shg;
  q[1] = sha * chb * chg - cha * shb * shg;
  q[2] = cha * shb * chg + sha * chb * shg;
  q[3] = cha * chb * shg - sha * shb * chg;

  itk::Versor<double> v;
  v.Set(q[0], q[1], q[2], q[3]);
  return v;
}

VersorTransformType::Pointer ConvertToVersorRigid3D(RigidTransformType::Pointer RT)
{
  VersorTransformType::Pointer VT = VersorTransformType::New();

  VT->SetFixedParameters( RT->GetFixedParameters() );

  itk::Matrix<double, 3, 3>           R = RT->GetMatrix();
  RigidTransformType::TranslationType T = RT->GetTranslation();

  VersorTransformType::ParametersType p;
  p.SetSize(6);
  itk::Versor<double> v;
  v.Set(R);
  // Get the first 3 elements of the versor;
  p[0] = v.GetRight()[0];
  p[1] = v.GetRight()[1];
  p[2] = v.GetRight()[2];
  p[3] = T[0];
  p[4] = T[1];
  p[5] = T[2];
  VT->SetParameters(p);
  return VT;
}

itk::Matrix<double, 3, 3> GetMatrixInverse(const itk::Matrix<double, 3, 3> & input)
{
  // Hack until enhancemetn fix in ITK is made public for assignment of matrix
  itk::Matrix<double, 3, 3>::InternalMatrixType temp = input.GetInverse();
  itk::Matrix<double, 3, 3>                     output;

  for( unsigned int r = 0; r < 3; ++r )
    {
    for( unsigned int c = 0; c < 3; ++c )
      {
      output(r, c) = temp(r, c);
      }
    }
  return output;
}

// Need to use a Pyramid filter here instead of just downsampling to a 32x32
// image
typedef itk::MultiResolutionPyramidImageFilter<SImageType, SImageType> PyramidFilterType;

PyramidFilterType::Pointer MakeThreeLevelPyramid(SImageType::Pointer refImage)
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
  for( unsigned int c = 0; c < pyramidSchedule.cols(); ++c )
    {
    // about 8mm
    pyramidSchedule[0][c] = static_cast<unsigned int>( 2 * round(4.0 / refImageSpacing[c]) );
    // about 4mm
    pyramidSchedule[1][c] = static_cast<unsigned int>( 2 * round(2.0 / refImageSpacing[c]) );
    // about 2mm
    pyramidSchedule[2][c] = static_cast<unsigned int>( 2 * round(1.0 / refImageSpacing[c]) );
    }
  MyPyramid->SetSchedule(pyramidSchedule);
  MyPyramid->Update();
  return MyPyramid;
}

PyramidFilterType::Pointer MakeOneLevelPyramid(SImageType::Pointer refImage)
{
  PyramidFilterType::ScheduleType pyramidSchedule;

  PyramidFilterType::Pointer MyPyramid = PyramidFilterType::New();

  MyPyramid->SetInput(refImage);
  MyPyramid->SetNumberOfLevels(1);
  pyramidSchedule.SetSize(1, 3);
  // Attempt to set a schedule so that the top of the pyramid
  // has images of about 8mm, and the next level has resolutions about 4mm
  // isotropic voxels
  // these are sizes found to work well for estimating MSP without making the
  // image too small.
  SImageType::SpacingType refImageSpacing = refImage->GetSpacing();
  for( unsigned int c = 0; c < pyramidSchedule.cols(); ++c )
    {
    pyramidSchedule[0][c] = static_cast<unsigned int>( vcl_floor(4.0 / refImageSpacing[c] + 0.5) );
    }
  MyPyramid->SetSchedule(pyramidSchedule);
  MyPyramid->Update();
  return MyPyramid;
}

// ////////////////////////////////////////////////////////////////////////
//  This is a lightweight wrapper for FindCenterOfBrainBasedOnTopOfHead, which
// evolved from TrimForegroundInDirection.
SImageType::PointType GetCenterOfHeadMass(SImageType::Pointer volume)
{
  // NOTE:   This is a bit of a hack to deal with vastly different sized heads,
  // note:   or where the field of view is not centered in the brain
  // note:   or where there is a massive amount of neck.

  // Choose the Inferior/Superior based on maximum dirction cosign

  // Find center of the image space which is where the MSP has been placed.
  const unsigned int    ISdirectionIndex = 2;
  SImageType::PointType CenterOfMass =
    FindCenterOfBrainBasedOnTopOfHead(volume,
                                      ISdirectionIndex,
                                      0.01, 7,
                                      700, 0);

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
RigidTransformType::Pointer computeTmspFromPoints(
  SImageType::PointType RP,
  SImageType::PointType AC,
  SImageType::PointType PC,
  SImageType::PointType DesiredCenter)
{
  // a variable to store correlation coefficient values
  SImageType::PointType::VectorType ACPC = PC - AC;

  ACPC.Normalize();
  SImageType::PointType::VectorType ACRP = RP - AC;

  vnl_vector_fixed<double, 3> NormalToMSPPlane = vnl_cross_3d( ACPC.GetVnlVector(), ACRP.GetVnlVector() );
  // --std::cout << ACPC << "=" << PC << " - " << AC << std::endl;
  // --std::cout << ACRP << "=" << RP << " - " << AC << std::endl;
  // --std::cout << NormalToMSPPlane << std::endl;
  SImageType::PointType::VectorType NormalToMSPPlaneVector;
  NormalToMSPPlaneVector[0] = NormalToMSPPlane[0];
  NormalToMSPPlaneVector[1] = NormalToMSPPlane[1];
  NormalToMSPPlaneVector[2] = NormalToMSPPlane[2];
  NormalToMSPPlaneVector.Normalize();
  if( NormalToMSPPlaneVector[0] < 0 )
    {
    NormalToMSPPlaneVector *= -1.0;
    }

  // --std::cout << "Norm: " << NormalToMSPPlaneVector << "     NORMACPC: " <<
  // ACPC << std::endl;

  // double
  // PlaneNormalBank=-vcl_atan2(NormalToMSPPlaneVector[2],NormalToMSPPlaneVector[0]);
  //  //Rotate the "Y" (i.e. Anterior to Posterior Axis)
  // double PlaneNormalHeading=-vcl_acos(NormalToMSPPlaneVector[0]);
  //               //Rotate the "Z" (i.e. Inferior to Superior Axis)

  double PlaneNormalBank = -vcl_atan2(NormalToMSPPlaneVector[2], NormalToMSPPlaneVector[0]);
  // Rotate the "Y" (i.e. Anterior to Posterior Axis)
  double PlaneNormalHeading = vcl_sin(NormalToMSPPlaneVector[1]);
  // Rotate the "Z" (i.e. Inferior to Superior Axis)
  double PlaneNormalAttitude = vcl_sin(ACPC[2]);

  SImageType::PointType::VectorType CenterOffset =
    AC.GetVectorFromOrigin()
    - DesiredCenter.GetVectorFromOrigin();

  RigidTransformType::Pointer AlignMSPTransform = RigidTransformType::New();
  AlignMSPTransform->SetCenter(DesiredCenter);
  AlignMSPTransform->SetRotation(PlaneNormalAttitude, PlaneNormalBank, PlaneNormalHeading);
    {
    // Clean up the rotation to make it orthogonal:
    const itk::Matrix<double, 3, 3> & CleanedOrthogonalized = itk::Orthogonalize3DRotationMatrix(
        AlignMSPTransform->GetMatrix() );
    AlignMSPTransform->SetMatrix( CleanedOrthogonalized );
    }
  AlignMSPTransform->SetTranslation(CenterOffset);
  if( LMC::globalverboseFlag )
    {
    std::cout << "\n============================" << CenterOffset << " \n" << AlignMSPTransform << std::endl;
    std::cout << " \n" << AlignMSPTransform->GetTranslation() << std::endl;
    std::cout << "============================\n" << std::endl;
    }
  return AlignMSPTransform;
}

// //
//
//
// ////////////////////////////////////////////////////////////////////////////////////////////////
int computeTemplateSize(const int r, const int h)
{
  const double h_2 = h / 2;
  const double r2 = r * r;
  int          size = 0;

  for( double k = -r; k <= r; ++k )
    {
    for( double j = -r; j <= r; ++j )
      {
      for( double i = -h_2; i <= h_2; ++i )
        {
        if( ( i * i + j * j )  <= r2 )
          {
          ++size;
          }
        }
      }
    }
  return size;
}

//
//
// ////////////////////////////////////////////////////////////////////////////////////////////////

// Use the law of cosines to determine the angle between the RPtoPC vector and
// the RPtoAC vector.
void decomposeRPAC(const SImageType::PointType & RP,
                   const SImageType::PointType & PC,
                   const SImageType::PointType & AC,
                   float *const RPPC_to_RPAC_angle, float *const RPAC_over_RPPC)
{
  const double RPtoPC = vcl_sqrt( ( ( RP[1] - PC[1] ) * ( RP[1] - PC[1] ) + ( RP[2] - PC[2] ) * ( RP[2] - PC[2] ) ) );
  const double RPtoAC = vcl_sqrt( ( ( RP[1] - AC[1] ) * ( RP[1] - AC[1] ) + ( RP[2] - AC[2] ) * ( RP[2] - AC[2] ) ) );
  const double PCtoAC = vcl_sqrt( ( ( PC[1] - AC[1] ) * ( PC[1] - AC[1] ) + ( PC[2] - AC[2] ) * ( PC[2] - AC[2] ) ) );

  const double cos_gamma = ( RPtoPC * RPtoPC + RPtoAC * RPtoAC - PCtoAC * PCtoAC ) / ( 2.0 * RPtoPC * RPtoAC );
  const double gamma = vcl_acos(cos_gamma);

  // const double sin_gamma = vcl_sqrt( 1.0 - cos_gamma*cos_gamma );

  *RPPC_to_RPAC_angle = gamma;       // The angle between the RP-PC vector and
                                     // the RP-AC vector
  *RPAC_over_RPPC = RPtoAC / RPtoPC; // The ratio of lengths between the RP-PC
                                     // vector and the RP-AC vector
                                     // std::cout << "----------TEST ANGLE " <<
                                     //  *RPPC_to_RPAC_angle *180.0/vnl_math::pi
                                     // << "  RATIO " << *RPAC_over_RPPC <<
                                     // std::endl;
}

SImageType::PointType::VectorType initialAC(const SImageType::PointType & RP,  const SImageType::PointType & PC,
                                            const double RPPC_to_RPAC_angleMean, const double RPAC_over_RPPCMean)
{
  const SImageType::PointType::VectorType RPPCVector = PC - RP;
  SImageType::PointType::VectorType       GuessAC;

  GuessAC[0] = ( RP[0] + RPPCVector[0] * .05 );
  const double cos_gamma = vcl_cos(RPPC_to_RPAC_angleMean);
  const double sin_gamma = vcl_sin(RPPC_to_RPAC_angleMean);
  // First rotate the RPPC vector in the direction of the AC poinnt
  GuessAC[1] = RP[1] + ( cos_gamma * RPPCVector[1] - sin_gamma * RPPCVector[2] ) * RPAC_over_RPPCMean;
  GuessAC[2] = RP[2] + ( sin_gamma * RPPCVector[1] + cos_gamma * RPPCVector[2] ) * RPAC_over_RPPCMean;
  return GuessAC;
}

void
extractArray(SImageType::Pointer & image,
             const SImageType::PointType & CenterPoint,
             const landmarksConstellationModelIO::IndexLocationVectorType & model,
             std::vector<float> & result_array)
{
  // TODO:  Move the LinearInterpolator to the function parameter, because the
  // image is never needed.
  LinearInterpolatorType::Pointer imInterp = LinearInterpolatorType::New();

  imInterp->SetInputImage(image);
  int q = 0;
  result_array.resize( model.size() );
  for( landmarksConstellationModelIO::IndexLocationVectorType::const_iterator it = model.begin();
       it != model.end(); ++it, ++q )
    {
    const SImageType::PointType point = CenterPoint + *it;
    bool                        isin = imInterp->IsInsideBuffer(point);
    if( isin )
      {
      result_array[q] = imInterp->Evaluate(point);
      }
    else
      {
      result_array[q] = 0.0;
      }
    }
}

void
extractArrayRemoveVectorMeanNormalize
  (SImageType::Pointer & image,
  const SImageType::PointType & CenterPoint,
  const landmarksConstellationModelIO::IndexLocationVectorType & model,
  std::vector<float> & result_array)
{
  // TODO:  Move the LinearInterpolator to the function parameter, because the
  // image is never needed.
  LinearInterpolatorType::Pointer imInterp = LinearInterpolatorType::New();

  imInterp->SetInputImage(image);
  result_array.resize( model.size() );
  double mean(0.0);
  landmarksConstellationModelIO::
  IndexLocationVectorType::const_iterator it = model.begin();
  std::vector<float>::iterator            resultIt = result_array.begin();
  // interpolate, accumulating mean as you go
  for( ; it != model.end(); ++it, ++resultIt )
    {
    const SImageType::PointType point = CenterPoint + *it;
    bool                        isin = imInterp->IsInsideBuffer(point);
    if( isin )
      {
      ( *resultIt ) = imInterp->Evaluate(point);
      }
    else
      {
      ( *resultIt ) = 0.0;
      }
    mean += ( *resultIt );
    }

  // remove mean, accumulating norm as you go.
  mean /= result_array.size();
  double norm = 0.0;
  for( resultIt = result_array.begin();
       resultIt != result_array.end(); ++resultIt )
    {
    double v( ( *resultIt ) - mean );
    *resultIt = v;
    norm += v * v;

    // debug abnormal behavior of norm!
    // std::cout << v << " " << norm << std::endl;
    }
  norm = vcl_sqrt(norm);
  if( norm < vcl_numeric_limits<double>::epsilon() )
    {
    std::cout << "WARNING:  ZERO NORM VECTOR." << __FILE__ << __LINE__ << std::endl;
    return;
    }
  // normalize array
  for( resultIt = result_array.begin();
       resultIt != result_array.end(); ++resultIt )
    {
    *resultIt /= norm;
    }
}

unsigned char
ShortToUChar(short in, short min, short max)
{
  double x(in);

  x -= static_cast<double>( min );
  double divisor = max - min;
  x /= divisor;
  x *= 255.0;
  return static_cast<unsigned char>( x );
}

// TODO:  Need to make these pointers const
SImageType::Pointer CreateTestCenteredRotatedImage2(const RigidTransformType::Pointer ACPC_MSP_AlignedTransform,
                                                    /* const
                                                      SImageType::PointType
                                                      finalPoint_MSP, */
                                                    const SImageType::PointType itkNotUsed( PreMSP_Point ),
                                                    /*const*/ SImageType::Pointer & image,
                                                    const RigidTransformType::Pointer & Point_Rotate)
{
  // ////// Compose test rotation with the translated ACPC alignment
  RigidTransformType::Pointer Point_Centered_TestRotated = RigidTransformType::New();

  Point_Centered_TestRotated->SetParameters( ACPC_MSP_AlignedTransform->GetParameters() );
  Point_Centered_TestRotated->SetFixedParameters( ACPC_MSP_AlignedTransform->GetFixedParameters() );
  Point_Centered_TestRotated->Compose(Point_Rotate, true); // TODO:  Perhaps
                                                           // change sign of
                                                           // Point_Translate,
                                                           // and the remove the
                                                           // true flag.

#if 1
  SImageType::Pointer image_Point_TestRotated = TransformResample<SImageType, SImageType>(
      image, image, /* DEFAULT TO ZERO */ 0,
      GetInterpolatorFromString<SImageType>("Linear"), Point_Centered_TestRotated.GetPointer() );
#else
  SImageType::Pointer image_Point_TestRotated = GenericImageResample<SImageType, SImageType, LinearInterpolatorType>(
      image,
      Point_Centered_TestRotated,
      /* DEFAULT TO ZERO */ 0);
#endif

  RigidTransformType::Pointer invPoint_Centered_TestRotated = RigidTransformType::New();
  Point_Centered_TestRotated->GetInverse(invPoint_Centered_TestRotated);
  return image_Point_TestRotated;
}

// typedef landmarksDataSet::PointType PointType;

void MakeLabelImage(SImageType::Pointer in,
                    const SImageType::PointType & RP,
                    const SImageType::PointType & AC,
                    const SImageType::PointType & PC,
                    const SImageType::PointType & VN4,
                    const std::string & fname)
{
  SImageType::Pointer maskImage = SImageType::New();

  maskImage->CopyInformation(in);
  maskImage->SetRegions( in->GetLargestPossibleRegion() );
  maskImage->Allocate();
  maskImage->FillBuffer(0);
    {
    SImageType::IndexType Index;
    maskImage->TransformPhysicalPointToIndex(RP, Index);
    maskImage->SetPixel(Index,  250);
    maskImage->TransformPhysicalPointToIndex(AC, Index);
    maskImage->SetPixel(Index,  255);
    maskImage->TransformPhysicalPointToIndex(PC, Index);
    maskImage->SetPixel(Index,  245);
    maskImage->TransformPhysicalPointToIndex(VN4, Index);
    maskImage->SetPixel(Index,  240);
    }
  itkUtil::WriteImage<SImageType>(maskImage, fname);
}

SImageType::PointType GetImageCenterPhysicalPoint(SImageType::Pointer & image)
{
  const SImageType::SizeType imageOverallSize = image->GetLargestPossibleRegion().GetSize();

  itk::ContinuousIndex<double, 3> centerIndex;

  for( int q = 0; q < SImageType::ImageDimension; ++q )
    {
    centerIndex[q] = 0.5 * ( imageOverallSize[q] - 1 );
    }
  SImageType::PointType centerLocation;
  image->TransformContinuousIndexToPhysicalPoint(centerIndex,  centerLocation);
  return centerLocation;
}

// TODO:  Need to make size and spacing command line arguments that are sent
// into
// the ResampleToIsotropicImage filter.
// HACK:  Currently hard-coded to 1.0^3 isotropic voxels in a 256^3 matrix
// --outputVolumeSize 256,256,256 --outputVolumeSpacing 1.0,1.0,1.0
SImageType::Pointer MakeIsoTropicReferenceImage(void)
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
