/*
*/
#ifndef __BRAINSConstellationDetectorPrimary__h
#define __BRAINSConstellationDetectorPrimary__h

#include "itksys/SystemTools.hxx"
#include "BRAINSThreadControl.h"

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

#include "itkResampleInPlaceImageFilter.h"
#include "itkMultiplyImageFilter.h"
#include "itkLandmarkBasedTransformInitializer.h"
#include "itkCastImageFilter.h"
#include <BRAINSFitHelper.h>

#include <cstdio>
#include <iostream>
#include <fstream>
#include <cstring>
#include <map>

class BRAINSConstellationDetectorPrimary
{
public:

  // Image, filter, transform typedef
  typedef short                                 PixelType;
  typedef itk::Image<PixelType, 3>              ImageType;
  typedef ImageType::Pointer                    ImagePointerType;
  typedef ImageType::PointType                  ImagePointType;
  typedef ImageType::SpacingType                ImageSpacingType;
  typedef ImageType::SizeType                   ImageSizeType;
  typedef ImageType::DirectionType              ImageDirectionType;
  typedef ImageType::IndexType                  ImageIndexType;
  typedef std::map<std::string, ImagePointType> LandmarksMapType;
  typedef std::map<std::string, double>         LandmarksWeightMapType;        // SHOULD BE DELETED

  typedef itk::ImageFileReader<ImageType>                         ReaderType;
  typedef itk::ImageFileWriter<ImageType>                         WriterType;
  typedef itk::FindCenterOfBrainFilter<ImageType>                 FindCenterFilter;
  typedef itk::BRAINSHoughEyeDetector<ImageType, ImageType>       HoughEyeDetectorType;
  typedef itk::BRAINSConstellationDetector2<ImageType, ImageType> Constellation2Type;
  typedef itk::TransformFileWriter                                TransformWriterType;
  typedef itk::VersorRigid3DTransform<double>                     VersorTransformType;
  std::string pathOut;
  std::string errorMsg;

  BRAINSConstellationDetectorPrimary();

  void SetHoughEyeDetectorMode(int houghEyeDetectorMode)
  {
    this->m_houghEyeDetectorMode = houghEyeDetectorMode;
  }

  void SetMspQualityLevel(unsigned int mspQualityLevel)
  {
    this->m_mspQualityLevel = mspQualityLevel;
  }

  void SetWritedebuggingImagesLevel(unsigned int writedebuggingImagesLevel)
  {
    this->m_writedebuggingImagesLevel = writedebuggingImagesLevel;
  }

  void SetNumberOfThreads(unsigned int numberOfThreads)
  {
    this->m_numberOfThreads = numberOfThreads;
  }

  void SetOtsuPercentileThreshold(double otsuPercentileThreshold)
  {
    this->m_otsuPercentileThreshold = otsuPercentileThreshold;
  }

  void SetAcLowerBound(double acLowerBound)
  {
    this->m_acLowerBound = acLowerBound;
  }

  void SetTrimRescaledIntensities(double trimRescaledIntensities)
  {
    this->m_trimRescaledIntensities = trimRescaledIntensities;
  }

  void SetRadiusMPJ(double radiusMPJ)
  {
    this->m_radiusMPJ = radiusMPJ;
  }

  void SetRadiusAC(double radiusAC)
  {
    this->m_radiusAC = radiusAC;
  }

  void SetRadiusPC(double radiusPC)
  {
    this->m_radiusPC = radiusPC;
  }

  void SetRadiusVN4(double radiusVN4)
  {
    this->m_radiusVN4 = radiusVN4;
  }

  void SetCutOutHeadInOutputVolume(bool cutOutHeadInOutputVolume)
  {
    this->m_cutOutHeadInOutputVolume = cutOutHeadInOutputVolume;
  }

  void SetRescaleIntensities(bool rescaleIntensities)
  {
    this->m_rescaleIntensities = rescaleIntensities;
  }

  void SetForceHoughEyeDetectorReportFailure(bool forceHoughEyeDetectorReportFailure)
  {
    this->m_forceHoughEyeDetectorReportFailure = forceHoughEyeDetectorReportFailure;
  }

  void SetDebug(bool debug)
  {
    this->m_debug = debug;
  }

  void SetVerbose(bool verbose)
  {
    this->m_verbose = verbose;
  }

  void SetAtlasVolume( const std::string & atlasVolume )
  {
    this->m_atlasVolume = atlasVolume;
  }

  void SetAtlasLandmarks ( const std::string & atlasLandmarks )
  {
    this->m_atlasLandmarks = atlasLandmarks;
  }

  void SetAtlasLandmarkWeights ( const std::string & atlasLandmarkWeights )
  {
    this->m_atlasLandmarkWeights = atlasLandmarkWeights;
  }

  void SetInputTemplateModel(std::string inputTemplateModel)
  {
    this->m_inputTemplateModel = inputTemplateModel;
  }

  void SetLLSModel(std::string llsModel)
  {
    this->m_llsModel = llsModel;
  }

  void SetInputVolume(std::string inputVolume)
  {
    this->m_inputVolume = inputVolume;
  }

  void SetOutputVolume(std::string outputVolume)
  {
    this->m_outputVolume = outputVolume;
  }

  void SetOutputResampledVolume(std::string outputResampledVolume)
  {
    this->m_outputResampledVolume = outputResampledVolume;
  }

  void SetOutputTransform(std::string outputTransform)
  {
    this->m_outputTransform = outputTransform;
  }

  void SetOutputLandmarksInInputSpace(std::string outputLandmarksInInputSpace)
  {
    this->m_outputLandmarksInInputSpace = outputLandmarksInInputSpace;
  }

  void SetOutputLandmarksInACPCAlignedSpace(std::string outputLandmarksInACPCAlignedSpace)
  {
    this->m_outputLandmarksInACPCAlignedSpace = outputLandmarksInACPCAlignedSpace;
  }

  void SetOutputLandmarkWeights(std::string outputLandmarkWeights)   // SHOULD BE DELETED
  {
    this->m_outputLandmarkWeights = outputLandmarkWeights;
  }

  void SetInputLandmarksPaired(std::string inputLandmarksPaired)     // SHOULD BE DELETED
  {
    this->m_inputLandmarksPaired = inputLandmarksPaired;
  }

  void SetOutputLandmarksPaired(std::string outputLandmarksPaired)   // SHOULD BE DELETED
  {
    this->m_outputLandmarksPaired = outputLandmarksPaired;
  }

  void SetOutputMRML(std::string outputMRML)
  {
    this->m_outputMRML = outputMRML;
  }

  void SetOutputVerificationScript(std::string outputVerificationScript)
  {
    this->m_outputVerificationScript = outputVerificationScript;
  }

  void SetOutputUntransformedClippedVolume(std::string outputUntransformedClippedVolume)
  {
    this->m_outputUntransformedClippedVolume = outputUntransformedClippedVolume;
  }

  void SetInputLandmarksEMSP(std::string inputLandmarksEMSP)
  {
    this->m_inputLandmarksEMSP = inputLandmarksEMSP;
  }

  void SetWriteBranded2DImage(std::string writeBranded2DImage)
  {
    this->m_writeBranded2DImage = writeBranded2DImage;
  }

  void SetBackgroundFillValueString(std::string backgroundFillValueString)
  {
    this->m_backgroundFillValueString = backgroundFillValueString;
  }

  void SetInterpolationMode(std::string interpolationMode)
  {
    this->m_interpolationMode = interpolationMode;
  }

  void SetRescaleIntensitiesOutputRange(std::vector<int> rescaleIntensitiesOutputRange)
  {
    this->m_rescaleIntensitiesOutputRange = rescaleIntensitiesOutputRange;
  }

  void SetForceACPoint(std::vector<float> forceACPoint)
  {
    this->m_forceACPoint = forceACPoint;
  }

  void SetForcePCPoint(std::vector<float> forcePCPoint)
  {
    this->m_forcePCPoint = forcePCPoint;
  }

  void SetForceVN4Point(std::vector<float> forceVN4Point)
  {
    this->m_forceVN4Point = forceVN4Point;
  }

  void SetForceRPPoint(std::vector<float> forceRPPoint)
  {
    this->m_forceRPPoint = forceRPPoint;
  }

  void SetResultsDir(std::string resultsDir)
  {
    this->m_resultsDir = resultsDir;
  }

  const LandmarksMapType & GetOutputLandmarksInInputSpace(void)
  {
    return this->m_outputLandmarksInInputSpaceMap;
  }

  const LandmarksMapType & GetOutputLandmarksInACPCAlignedSpace(void)
  {
    return this->m_outputLandmarksInACPCAlignedSpaceMap;
  }

  bool Compute(void);

private:

  int          m_houghEyeDetectorMode;      // 1
  unsigned int m_mspQualityLevel;           // 2
  unsigned int m_writedebuggingImagesLevel; // 0
  unsigned int m_numberOfThreads;           // 1
  double       m_otsuPercentileThreshold;   // 0.01
  double       m_acLowerBound;              // 1000.0
  double       m_trimRescaledIntensities;   // 4.4172

  double m_radiusMPJ;   // -1
  double m_radiusAC;    // -1
  double m_radiusPC;    // -1
  double m_radiusVN4;   // -1

  bool m_cutOutHeadInOutputVolume;              // false
  bool m_rescaleIntensities;                    // false
  bool m_forceHoughEyeDetectorReportFailure;    // false
  bool m_debug;                                 // false
  bool m_verbose;                               // false

  std::string m_inputTemplateModel;
  std::string m_llsModel;
  std::string m_inputVolume;
  std::string m_outputVolume;
  std::string m_outputResampledVolume;
  std::string m_outputTransform;
  std::string m_outputLandmarksInInputSpace;
  std::string m_outputLandmarksInACPCAlignedSpace;
  std::string m_outputLandmarkWeights;   // SHOULD BE DELETED
  std::string m_inputLandmarksPaired;    // SHOULD BE DELETED
  std::string m_outputLandmarksPaired;   // SHOULD BE DELETED
  std::string m_outputMRML;
  std::string m_outputVerificationScript;
  std::string m_outputUntransformedClippedVolume;
  std::string m_inputLandmarksEMSP;
  std::string m_writeBranded2DImage;
  std::string m_backgroundFillValueString;
  std::string m_interpolationMode;
  std::string m_atlasVolume ;
  std::string m_atlasLandmarks ;
  std::string m_atlasLandmarkWeights ;

  std::vector<int> m_rescaleIntensitiesOutputRange;   // default = [40,4000]

  std::vector<float> m_forceACPoint;    // default = 0.
  std::vector<float> m_forcePCPoint;    // default = 0.
  std::vector<float> m_forceVN4Point;   // default = 0.
  std::vector<float> m_forceRPPoint;    // default = 0.

  std::string m_resultsDir;   // default = "./"

  LandmarksMapType m_outputLandmarksInInputSpaceMap;
  LandmarksMapType m_outputLandmarksInACPCAlignedSpaceMap;
};

#endif
