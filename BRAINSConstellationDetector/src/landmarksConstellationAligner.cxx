/*
 * Author: Ali Ghayoor
 * at Psychiatry Imaging Lab,
 * University of Iowa Health Care 2011
 */

// I N C L U D E S ////////////////////////////////////////////////////////////

#include "itksys/SystemTools.hxx"
#include "BRAINSThreadControl.h"
#include "Slicer3LandmarkIO.h"

#include "itkCommand.h"
#include "itkImage.h"
#include "itkOrthogonalize3DRotationMatrix.h"
#include "itkVersorRigid3DTransform.h"
#include "itkEuler3DTransform.h"

#include "stdio.h"
#include <iostream>
#include <fstream>
#include "string.h"
#include <map>
#include <vnl/vnl_cross.h>

#include "landmarksConstellationAlignerCLP.h"

// D E F I N E S //////////////////////////////////////////////////////////////

const unsigned int LocalImageDimension = 3;
typedef short PixelType;

typedef itk::Image<PixelType, LocalImageDimension> SImageType;
typedef SImageType::PointType                      ImagePointType;
typedef std::map<std::string, ImagePointType>      LandmarksMapType;
typedef itk::VersorRigid3DTransform<double>        VersorTransformType;
typedef itk::Euler3DTransform<double>              RigidTransformType;

// F U N C T I O N S //////////////////////////////////////////////////////////

SImageType::PointType GetOriginalSpaceNamedPoint(const LandmarksMapType & landmarks, const std::string & NamedPoint)
{
  LandmarksMapType::const_iterator itpair = landmarks.find(NamedPoint);

  if( itpair == landmarks.end() )
    {
    std::cout << "ERROR:  " << NamedPoint << " not found in list." << std::endl;
    return SImageType::PointType();
    }
  return itpair->second;
}

//////////////////
RigidTransformType::Pointer computeTmspFromPoints(SImageType::PointType RP,
                                                  SImageType::PointType AC,
                                                  SImageType::PointType PC,
                                                  SImageType::PointType DesiredCenter)
{
  // a variable to store correlation coefficient values
  SImageType::PointType::VectorType ACPC = PC - AC;

  ACPC.Normalize();
  SImageType::PointType::VectorType ACRP = RP - AC;

  vnl_vector_fixed<double, 3>       NormalToMSPPlane = vnl_cross_3d( ACPC.GetVnlVector(), ACRP.GetVnlVector() );
  SImageType::PointType::VectorType NormalToMSPPlaneVector;
  NormalToMSPPlaneVector[0] = NormalToMSPPlane[0];
  NormalToMSPPlaneVector[1] = NormalToMSPPlane[1];
  NormalToMSPPlaneVector[2] = NormalToMSPPlane[2];
  NormalToMSPPlaneVector.Normalize();
  if( NormalToMSPPlaneVector[0] < 0 )
    {
    NormalToMSPPlaneVector *= -1.0;
    }

  double PlaneNormalBank = -vcl_atan2(NormalToMSPPlaneVector[2], NormalToMSPPlaneVector[0]);
  // Rotate the "Y" (i.e. Anterior to Posterior Axis)
  double PlaneNormalHeading = vcl_sin(NormalToMSPPlaneVector[1]);
  // Rotate the "Z" (i.e. Inferior to Superior Axis)
  double PlaneNormalAttitude = vcl_sin(ACPC[2]);

  SImageType::PointType::VectorType CenterOffset = AC.GetVectorFromOrigin() - DesiredCenter.GetVectorFromOrigin();

  RigidTransformType::Pointer AlignMSPTransform = RigidTransformType::New();
  AlignMSPTransform->SetCenter(DesiredCenter);
  AlignMSPTransform->SetRotation(PlaneNormalAttitude, PlaneNormalBank, PlaneNormalHeading);
    {
    // Clean up the rotation to make it orthogonal:
    const itk::Matrix<double, 3, 3> & CleanedOrthogonalized = itk::Orthogonalize3DRotationMatrix(
        AlignMSPTransform->GetRotationMatrix() );
    AlignMSPTransform->SetRotationMatrix( CleanedOrthogonalized );
    }
  AlignMSPTransform->SetTranslation(CenterOffset);

  return AlignMSPTransform;
}

//////////////////
RigidTransformType::Pointer GetACPCAlignedZeroCenteredTransform(const LandmarksMapType & landmarks)
{
  SImageType::PointType ZeroCenter;

  ZeroCenter.Fill(0.0);
  RigidTransformType::Pointer landmarkDefinedACPCAlignedToZeroTransform = computeTmspFromPoints(GetOriginalSpaceNamedPoint(
                                                                                                  landmarks,
                                                                                                  "RP"),
                                                                                                GetOriginalSpaceNamedPoint(
                                                                                                  landmarks,
                                                                                                  "AC"),
                                                                                                GetOriginalSpaceNamedPoint(
                                                                                                  landmarks, "PC"),
                                                                                                ZeroCenter);
  return landmarkDefinedACPCAlignedToZeroTransform;
}

// M A I N ////////////////////////////////////////////////////////////////////

int main( int argc, char *argv[] )
{
  PARSE_ARGS;

  // Verify input parameters
  if( ( ( outputLandmarksPaired.compare( "" ) != 0 ) &&
        ( inputLandmarksPaired.compare( "" ) == 0 ) )
      ||
      ( ( outputLandmarksPaired.compare( "" ) == 0 ) &&
        ( inputLandmarksPaired.compare( "" ) != 0 ) ) )
    {
    std::cerr << "The outputLandmark parameter should be paired with"
              << "the inputLandmark parameter." << std::endl;
    std::cerr << "No output acpc-aligned landmark list file is generated" << std::endl;
    std::cerr << "Type " << argv[0] << " -h for more help." << std::endl;
    }

  if( ( outputLandmarksPaired.compare( "" ) != 0 ) &&
      ( inputLandmarksPaired.compare( "" ) != 0 ) )
    {
    LandmarksMapType origLandmarks = ReadSlicer3toITKLmk( inputLandmarksPaired );
    LandmarksMapType alignedLandmarks;

    // computing the forward and inverse transform
    RigidTransformType::Pointer ZeroCenteredTransform = GetACPCAlignedZeroCenteredTransform( origLandmarks );

    VersorTransformType::Pointer finalTransform = VersorTransformType::New();
    finalTransform->SetFixedParameters( ZeroCenteredTransform->GetFixedParameters() );
    itk::Versor<double>               versorRotation;
    const itk::Matrix<double, 3, 3> & CleanedOrthogonalized = itk::Orthogonalize3DRotationMatrix(
        ZeroCenteredTransform->GetRotationMatrix() );
    versorRotation.Set( CleanedOrthogonalized );
    finalTransform->SetRotation( versorRotation );
    finalTransform->SetTranslation( ZeroCenteredTransform->GetTranslation() );
    // inverse transform
    VersorTransformType::Pointer invFinalTransform = VersorTransformType::New();
    SImageType::PointType        centerPoint = finalTransform->GetCenter();
    invFinalTransform->SetCenter( centerPoint );
    invFinalTransform->SetIdentity();
    finalTransform->GetInverse( invFinalTransform );

    // converting the original landmark file to the aligned landmark file
    LandmarksMapType::const_iterator olit = origLandmarks.begin();
    for( ; olit != origLandmarks.end(); ++olit )
      {
      alignedLandmarks[olit->first] = invFinalTransform->TransformPoint( olit->second );
      }
    // writing the acpc-aligned landmark file
    WriteITKtoSlicer3Lmk( outputLandmarksPaired, alignedLandmarks );
    std::cout << "The acpc-aligned landmark list file is written." << std::endl;
    }

  return 0;
}
