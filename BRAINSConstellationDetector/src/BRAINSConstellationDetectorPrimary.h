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
 */
#ifndef __BRAINSConstellationDetectorPrimary__h
#define __BRAINSConstellationDetectorPrimary__h

#include "itksys/SystemTools.hxx"
#include "BRAINSThreadControl.h"

#include "StandardizeMaskIntensity.h"
#include "landmarksConstellationCommon.h"
#include "itkFindCenterOfBrainFilter.h"
#include "BRAINSHoughEyeDetector.h"
#include "BRAINSConstellationDetector2.h"
#include "Slicer3LandmarkIO.h"
#include "Slicer3LandmarkWeightIO.h"
#include "LLSModel.h"

#include "itkCommand.h"
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkTransformFileWriter.h"
#include "itkVersorRigid3DTransform.h"

#include <cstdio>
#include <iostream>
#include <fstream>
#include <cstring>
#include <map>

class BRAINSConstellationDetectorPrimary
{
public:
  // Image, filter, transform typedef
  using PixelType = short;
  using ImageType = itk::Image<PixelType, 3>;
  using ImagePointerType = ImageType::Pointer;
  using ImagePointType = ImageType::PointType;
  using ImageSpacingType = ImageType::SpacingType;
  using ImageSizeType = ImageType::SizeType;
  using ImageDirectionType = ImageType::DirectionType;
  using ImageIndexType = ImageType::IndexType;

  using ReaderType = itk::ImageFileReader<ImageType>;
  using WriterType = itk::ImageFileWriter<ImageType>;
  using FindCenterFilter = itk::FindCenterOfBrainFilter<ImageType>;
  using HoughEyeDetectorType = itk::BRAINSHoughEyeDetector<ImageType, ImageType>;
  using TransformWriterType = itk::TransformFileWriter;
  using VersorRigidTransformType = itk::VersorRigid3DTransform<double>;
  std::string pathOut;
  std::string errorMsg;

  BRAINSConstellationDetectorPrimary();

  void
  SetHoughEyeDetectorMode(int houghEyeDetectorMode)
  {
    this->m_houghEyeDetectorMode = houghEyeDetectorMode;
  }

  void
  SetMspQualityLevel(unsigned int mspQualityLevel)
  {
    this->m_mspQualityLevel = mspQualityLevel;
  }

  void
  SetWritedebuggingImagesLevel(unsigned int writedebuggingImagesLevel)
  {
    this->m_writedebuggingImagesLevel = writedebuggingImagesLevel;
  }

  void
  SetNumberOfWorkUnits(unsigned int numberOfThreads)
  {
    this->m_numberOfThreads = numberOfThreads;
  }

  void
  SetOtsuPercentileThreshold(double otsuPercentileThreshold)
  {
    this->m_otsuPercentileThreshold = otsuPercentileThreshold;
  }

  void
  SetAcLowerBound(double acLowerBound)
  {
    this->m_acLowerBound = acLowerBound;
  }

  void
  SetTrimRescaledIntensities(double trimRescaledIntensities)
  {
    this->m_trimRescaledIntensities = trimRescaledIntensities;
  }

  void
  SetRadiusMPJ(double radiusMPJ)
  {
    this->m_radiusMPJ = radiusMPJ;
  }

  void
  SetRadiusAC(double radiusAC)
  {
    this->m_radiusAC = radiusAC;
  }

  void
  SetRadiusPC(double radiusPC)
  {
    this->m_radiusPC = radiusPC;
  }

  void
  SetRadiusVN4(double radiusVN4)
  {
    this->m_radiusVN4 = radiusVN4;
  }

  void
  SetCutOutHeadInOutputVolume(bool cutOutHeadInOutputVolume)
  {
    this->m_cutOutHeadInOutputVolume = cutOutHeadInOutputVolume;
  }

  void
  SetRescaleIntensities(bool rescaleIntensities)
  {
    this->m_rescaleIntensities = rescaleIntensities;
  }

  void
  SetForceHoughEyeDetectorReportFailure(bool forceHoughEyeDetectorReportFailure)
  {
    this->m_forceHoughEyeDetectorReportFailure = forceHoughEyeDetectorReportFailure;
  }

  void
  SetDebug(bool debug)
  {
    this->m_debug = debug;
  }

  void
  SetVerbose(bool verbose)
  {
    this->m_verbose = verbose;
  }

  void
  SetAtlasVolume(const std::string & atlasVolume)
  {
    this->m_atlasVolume = atlasVolume;
  }

  void
  SetAtlasLandmarks(const std::string & atlasLandmarks)
  {
    this->m_atlasLandmarks = atlasLandmarks;
  }

  void
  SetAtlasLandmarkWeights(const std::string & atlasLandmarkWeights)
  {
    this->m_atlasLandmarkWeights = atlasLandmarkWeights;
  }

  void
  SetInputTemplateModel(std::string inputTemplateModel)
  {
    this->m_inputTemplateModel = inputTemplateModel;
  }

  void
  SetLLSModel(std::string llsModel)
  {
    this->m_llsModel = llsModel;
  }

  void
  SetInputVolume(std::string inputVolume)
  {
    this->m_inputVolume = inputVolume;
  }

  void
  SetOutputVolume(std::string outputVolume)
  {
    this->m_outputVolume = outputVolume;
  }

  void
  SetOutputResampledVolume(std::string outputResampledVolume)
  {
    this->m_outputResampledVolume = outputResampledVolume;
  }

  void
  SetOutputTransform(std::string outputTransform)
  {
    this->m_outputTransform = outputTransform;
  }

  void
  SetOutputLandmarksInInputSpace(std::string outputLandmarksInInputSpace)
  {
    this->m_outputLandmarksInInputSpace = outputLandmarksInInputSpace;
  }

  void
  SetOutputLandmarksInACPCAlignedSpace(std::string outputLandmarksInACPCAlignedSpace)
  {
    this->m_outputLandmarksInACPCAlignedSpace = outputLandmarksInACPCAlignedSpace;
  }

  void
  SetOutputMRML(std::string outputMRML)
  {
    this->m_outputMRML = outputMRML;
  }

  void
  SetOutputVerificationScript(std::string outputVerificationScript)
  {
    this->m_outputVerificationScript = outputVerificationScript;
  }

  void
  SetOutputUntransformedClippedVolume(std::string outputUntransformedClippedVolume)
  {
    this->m_outputUntransformedClippedVolume = outputUntransformedClippedVolume;
  }

  void
  SetInputLandmarksEMSP(std::string inputLandmarksEMSP)
  {
    this->orig_lmks_filename = inputLandmarksEMSP;
  }

  void
  SetWriteBranded2DImage(std::string writeBranded2DImage)
  {
    this->m_writeBranded2DImage = writeBranded2DImage;
  }

  void
  SetBackgroundFillValueString(std::string backgroundFillValueString)
  {
    this->m_backgroundFillValueString = backgroundFillValueString;
  }

  void
  SetInterpolationMode(std::string interpolationMode)
  {
    this->m_interpolationMode = interpolationMode;
  }

  void
  SetRescaleIntensitiesOutputRange(std::vector<int> rescaleIntensitiesOutputRange)
  {
    this->m_rescaleIntensitiesOutputRange = rescaleIntensitiesOutputRange;
  }


  void
  SetForce_orig_lmk_ACPointRAS(std::vector<float> forceACPointRAS)
  {
    this->m_force_orig_lmk_ACPointLPS = RAS2LPS(forceACPointRAS, "AC");
  }

  void
  SetForce_orig_lmk_PCPointRAS(std::vector<float> forcePCPointRAS)
  {
    this->m_force_orig_lmk_PCPointLPS = RAS2LPS(forcePCPointRAS, "PC");
  }

  void
  SetForce_with_lmk_VN4PointRAS(std::vector<float> forceVN4PointRAS)
  {
    this->m_force_orig_lmk_VN4PointLPS = RAS2LPS(forceVN4PointRAS, "VN4");
  }

  void
  SetForce_with_lmk_RPPointRAS(std::vector<float> forceRPPointRAS)
  {
    this->m_force_orig_lmk_RPPointLPS = RAS2LPS(forceRPPointRAS, "RP");
  }

  void
  SetResultsDir(std::string resultsDir)
  {
    this->m_resultsDir = resultsDir;
  }

  const LandmarksMapType &
  GetOutputLandmarksInInputSpace()
  {
    return this->m_outputLandmarksInInputSpaceMap;
  }

  const LandmarksMapType &
  GetOutputLandmarksInACPCAlignedSpace()
  {
    return this->m_outputLandmarksInACPCAlignedSpaceMap;
  }

  bool
  Compute();

private:
  static std::vector<float> &
  RAS2LPS(std::vector<float> & RASlmk, const std::string & name)
  {
    if (!RASlmk.empty())
    {
      // If not empty, then only size 3 is allowed, and must negate first two elements.
      if (RASlmk.size() != 3)
      {
        std::cerr << "Forced landmark requires 3 dimensions, " << RASlmk.size() << " provided for " << name << ". ";
        for (const auto e : RASlmk)
        {
          std::cerr << e << ", ";
        }
        std::cerr << std::endl;
        RASlmk.clear();
        exit(-1);
      }
      RASlmk[0] *= -1.0f;
      RASlmk[1] *= -1.0f;
    }
    return RASlmk;
  }
  static ImagePointType
  localFindCenterHeadFunc(const ImageType::ConstPointer & img);

  int          m_houghEyeDetectorMode;      // 1
  unsigned int m_mspQualityLevel;           // 2
  unsigned int m_writedebuggingImagesLevel; // 0
  unsigned int m_numberOfThreads;           // 1
  double       m_otsuPercentileThreshold;   // 0.01
  double       m_acLowerBound;              // 1000.0
  double       m_trimRescaledIntensities;   // 4.4172

  double m_radiusMPJ; // -1
  double m_radiusAC;  // -1
  double m_radiusPC;  // -1
  double m_radiusVN4; // -1

  bool m_cutOutHeadInOutputVolume;           // false
  bool m_rescaleIntensities;                 // false
  bool m_forceHoughEyeDetectorReportFailure; // false
  bool m_debug;                              // false
  bool m_verbose;                            // false

  std::string m_inputTemplateModel;
  std::string m_llsModel;
  std::string m_inputVolume;
  std::string m_outputVolume;
  std::string m_outputResampledVolume;
  std::string m_outputTransform;
  std::string m_outputLandmarksInInputSpace;
  std::string m_outputLandmarksInACPCAlignedSpace;
  std::string m_outputMRML;
  std::string m_outputVerificationScript;
  std::string m_outputUntransformedClippedVolume;
  std::string orig_lmks_filename;
  std::string m_writeBranded2DImage;
  std::string m_backgroundFillValueString;
  std::string m_interpolationMode;
  std::string m_atlasVolume;
  std::string m_atlasLandmarks;
  std::string m_atlasLandmarkWeights;

  std::vector<int> m_rescaleIntensitiesOutputRange; // default = [40,4000]

  std::vector<float> m_force_orig_lmk_ACPointLPS;  // default = 0.
  std::vector<float> m_force_orig_lmk_PCPointLPS;  // default = 0.
  std::vector<float> m_force_orig_lmk_VN4PointLPS; // default = 0.
  std::vector<float> m_force_orig_lmk_RPPointLPS;  // default = 0.

  std::string m_resultsDir; // default = "./"

  LandmarksMapType m_outputLandmarksInInputSpaceMap;
  LandmarksMapType m_outputLandmarksInACPCAlignedSpaceMap;
};

#endif
