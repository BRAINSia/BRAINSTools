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
//
//  main.cpp
//  sphericalCoordinateGenerationTest
//
//  Created by Ali Ghayoor on 1/14/13.
//  Copyright (c) 2013 Ali Ghayoor. All rights reserved.
//

#include "BRAINSCutDataHandler.h"
#include "BRAINSCutConfiguration.h"
#include <itkIO.h>
#include "SphericalCoordinateGenerationCLP.h"

// Define the required functions

WorkingImageType::IndexType::IndexValueType
TruncatedHalf(const WorkingImageType::SizeType::SizeValueType & v)
{
  return static_cast<WorkingImageType::IndexType::IndexValueType>(static_cast<double>(v) * 0.5);
}

void
CreateNewFloatImageFromTemplate(WorkingImageType::Pointer &       PointerToOutputImage,
                                const WorkingImageType::Pointer & PreInitializedImage)
{
  WorkingImageType::RegionType region;

  PointerToOutputImage = WorkingImageType::New();
  region.SetSize(PreInitializedImage->GetLargestPossibleRegion().GetSize());
  region.SetIndex(PreInitializedImage->GetLargestPossibleRegion().GetIndex());
  PointerToOutputImage->SetLargestPossibleRegion(region);
  PointerToOutputImage->SetBufferedRegion(region);
  PointerToOutputImage->SetRequestedRegion(region);
  PointerToOutputImage->CopyInformation(PreInitializedImage);
  PointerToOutputImage->Allocate();
  PointerToOutputImage->FillBuffer(0.0);
  // NO LONGER NEEDED CHECK_CORONAL(PointerToOutputImage->GetDirection());
  PointerToOutputImage->SetDirection(PreInitializedImage->GetDirection());
  PointerToOutputImage->SetMetaDataDictionary(PreInitializedImage->GetMetaDataDictionary());
  itk::ImageRegionIterator<WorkingImageType> bbri(PointerToOutputImage,
                                                  PointerToOutputImage->GetLargestPossibleRegion());
  bbri.GoToBegin();
  while (!bbri.IsAtEnd())
  {
    // Zeroing voxel signal intensity values
    bbri.Set(itk::NumericTraits<WorkingImageType::PixelType>::ZeroValue());
    ++bbri;
  }
}

void
XYZToSpherical(const itk::Point<float, 3> & LocationWithOriginAtCenterOfImage,
               float &                      rhoValue,
               float &                      phiValue,
               float &                      thetaValue)
{
  /*Rho*/
#define _SQR(a) ((a) * (a))
  rhoValue = static_cast<float>(std::sqrt(_SQR(LocationWithOriginAtCenterOfImage[0]) +
                                          _SQR(LocationWithOriginAtCenterOfImage[1]) +
                                          _SQR(LocationWithOriginAtCenterOfImage[2])));
#undef _SQR
  /*Phi*/
  phiValue = 0.0F;
  if (LocationWithOriginAtCenterOfImage[0] < 0)
  {
    phiValue = std::atan2(-LocationWithOriginAtCenterOfImage[0], LocationWithOriginAtCenterOfImage[1]);
  }
  else
  {
    phiValue = std::atan2(LocationWithOriginAtCenterOfImage[0], LocationWithOriginAtCenterOfImage[1]);
  }
  /*Theta*/
  thetaValue = 0.0F;
  if (LocationWithOriginAtCenterOfImage[2] < 0)
  {
    thetaValue = std::atan2(-LocationWithOriginAtCenterOfImage[2], LocationWithOriginAtCenterOfImage[1]);
  }
  else
  {
    thetaValue = std::atan2(LocationWithOriginAtCenterOfImage[2], LocationWithOriginAtCenterOfImage[1]);
  }

  //  thetaValue = std::acos(LocationWithOriginAtCenterOfImage[2]/rhoValue);

  rhoValue = rhoValue / 128.0F; // The largest brain ever will always fit in a sphere
  // with radius of 128MM centered at the AC point
  phiValue = phiValue / (itk::Math::pi);
  thetaValue = thetaValue / (itk::Math::pi);
}

// main program

int
main(int argc, const char * argv[])
{
  PARSE_ARGS;
  BRAINSRegisterAlternateIO();

  // read the atlas image from input
  WorkingImageType::Pointer atlasImage = itkUtil::ReadImage<WorkingImageType>(inputAtlasImage);

  const WorkingImageType::SizeType atlasImageSize(atlasImage->GetLargestPossibleRegion().GetSize());

  WorkingImageType::IndexType centerOfAtlas;

  centerOfAtlas[0] = TruncatedHalf(atlasImageSize[0]);
  centerOfAtlas[1] = TruncatedHalf(atlasImageSize[1]);
  centerOfAtlas[2] = TruncatedHalf(atlasImageSize[2]);

  itk::Point<WorkingPixelType, DIMENSION> centerOfAtlasPhysicalSpace;
  atlasImage->TransformIndexToPhysicalPoint(centerOfAtlas, centerOfAtlasPhysicalSpace);

  WorkingImageType::Pointer rhoImage, phiImage, thetaImage;

  CreateNewFloatImageFromTemplate(rhoImage, atlasImage);
  CreateNewFloatImageFromTemplate(phiImage, atlasImage);
  CreateNewFloatImageFromTemplate(thetaImage, atlasImage);

  itk::ImageRegionIterator<WorkingImageType> it(atlasImage, atlasImage->GetLargestPossibleRegion());
  it.GoToBegin();

  itk::ImageRegionIterator<WorkingImageType> rhoit(rhoImage, rhoImage->GetLargestPossibleRegion());
  rhoit.GoToBegin();

  itk::ImageRegionIterator<WorkingImageType> phiit(phiImage, phiImage->GetLargestPossibleRegion());
  phiit.GoToBegin();

  itk::ImageRegionIterator<WorkingImageType> thetait(thetaImage, thetaImage->GetLargestPossibleRegion());
  thetait.GoToBegin();

  itk::Point<float, DIMENSION> currentLocationPhysicalSpace;
  itk::Point<float, DIMENSION> LocationWithRespectToCenterOfImageInMM;

  while (!it.IsAtEnd())
  {
    const WorkingImageType::IndexType CurrentIndex = it.GetIndex();
    atlasImage->TransformIndexToPhysicalPoint(CurrentIndex, currentLocationPhysicalSpace);
    for (unsigned i = 0; i < DIMENSION; i++)
    {
      LocationWithRespectToCenterOfImageInMM[i] = currentLocationPhysicalSpace[i] - centerOfAtlasPhysicalSpace[i];
    }

    {
      float rhoValue, phiValue, thetaValue;
      XYZToSpherical(LocationWithRespectToCenterOfImageInMM, rhoValue, phiValue, thetaValue);
      rhoit.Set(rhoValue);
      phiit.Set(phiValue);
      thetait.Set(thetaValue);
    }

    ++it;
    ++rhoit;
    ++phiit;
    ++thetait;
  }

  std::string rhoPath = outputPath + "/rho.nii.gz";
  std::string phiPath = outputPath + "/phi.nii.gz";
  std::string thetaPath = outputPath + "/theta.nii.gz";

  // Write the output images to the disk.
  itkUtil::WriteImage<WorkingImageType>(rhoImage, rhoPath);
  itkUtil::WriteImage<WorkingImageType>(phiImage, phiPath);
  itkUtil::WriteImage<WorkingImageType>(thetaImage, thetaPath);

  return EXIT_SUCCESS;
}
