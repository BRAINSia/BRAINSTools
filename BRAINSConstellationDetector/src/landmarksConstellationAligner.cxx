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

#include "PrepareOutputImages.h"

#include <cstdio>
#include <iostream>
#include <fstream>
#include <cstring>
#include <map>
#include <vnl/vnl_cross.h>

#include "landmarksConstellationAlignerCLP.h"

// D E F I N E S //////////////////////////////////////////////////////////////
using ImagePointType = SImageType::PointType;

// F U N C T I O N S //////////////////////////////////////////////////////////
RigidTransformType::Pointer
GetACPCAlignedZeroCenteredTransform(const LandmarksMapType & landmarks)
{
  SImageType::PointType ZeroCenter;

  ZeroCenter.Fill(0.0);
  RigidTransformType::Pointer landmarkDefinedACPCAlignedToZeroTransform =
    computeTmspFromPoints(GetNamedPointFromLandmarkList(landmarks, "RP"),
                          GetNamedPointFromLandmarkList(landmarks, "AC"),
                          GetNamedPointFromLandmarkList(landmarks, "PC"),
                          ZeroCenter);
  return landmarkDefinedACPCAlignedToZeroTransform;
}

// M A I N ////////////////////////////////////////////////////////////////////
int
main(int argc, char * argv[])
{
  PARSE_ARGS;
  BRAINSRegisterAlternateIO();

  // Verify input parameters
  if (((outputLandmarksPaired.compare("") != 0) && (inputLandmarksPaired.compare("") == 0)) ||
      ((outputLandmarksPaired.compare("") == 0) && (inputLandmarksPaired.compare("") != 0)))
  {
    std::cerr << "The outputLandmark parameter should be paired with"
              << "the inputLandmark parameter." << std::endl;
    std::cerr << "No output acpc-aligned landmark list file is generated" << std::endl;
    std::cerr << "Type " << argv[0] << " -h for more help." << std::endl;
  }

  if ((outputLandmarksPaired.compare("") != 0) && (inputLandmarksPaired.compare("") != 0))
  {
    LandmarksMapType origLandmarks = ReadSlicer3toITKLmk(inputLandmarksPaired);
    LandmarksMapType alignedLandmarks;

    // computing the forward and inverse transform
    RigidTransformType::Pointer ZeroCenteredTransform = GetACPCAlignedZeroCenteredTransform(origLandmarks);

    VersorTransformType::Pointer finalTransform = VersorTransformType::New();
    finalTransform->SetFixedParameters(ZeroCenteredTransform->GetFixedParameters());
    itk::Versor<double>               versorRotation;
    const itk::Matrix<double, 3, 3> & CleanedOrthogonalized =
      itk::Orthogonalize3DRotationMatrix(ZeroCenteredTransform->GetMatrix());
    versorRotation.Set(CleanedOrthogonalized);
    finalTransform->SetRotation(versorRotation);
    finalTransform->SetTranslation(ZeroCenteredTransform->GetTranslation());
    // inverse transform
    VersorTransformType::Pointer invFinalTransform = VersorTransformType::New();
    SImageType::PointType        centerPoint = finalTransform->GetCenter();
    invFinalTransform->SetCenter(centerPoint);
    invFinalTransform->SetIdentity();
    finalTransform->GetInverse(invFinalTransform);

    // converting the original landmark file to the aligned landmark file
    LandmarksMapType::const_iterator olit = origLandmarks.begin();
    for (; olit != origLandmarks.end(); ++olit)
    {
      alignedLandmarks[olit->first] = invFinalTransform->TransformPoint(olit->second);
    }
    // writing the acpc-aligned landmark file
    WriteITKtoSlicer3Lmk(outputLandmarksPaired, alignedLandmarks);
    std::cout << "The acpc-aligned landmark list file is written." << std::endl;
  }

  return 0;
}
