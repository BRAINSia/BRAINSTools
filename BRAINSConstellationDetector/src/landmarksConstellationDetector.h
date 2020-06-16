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
 * Author: Han J. Johnson, Wei Lu  at SINAPSE Lab
 * University of Iowa Health Care 2010
 */

#ifndef __landmarksConstellationDetector__h
#define __landmarksConstellationDetector__h
#include <map>

#include "landmarksConstellationCommon.h"
#include "landmarksConstellationModelIO.h"
#include "itkAffineTransform.h"
#include "itkOtsuThresholdImageFilter.h"
#include "Slicer3LandmarkIO.h"
#include "PrepareOutputImages.h"

#include "itkMaskedFFTNormalizedCorrelationImageFilter.h"
#include "itkMinimumMaximumImageCalculator.h"
#include "itkBinaryImageToLabelMapFilter.h"
#include "itkLabelMapToLabelImageFilter.h"
#include "itkLabelStatisticsImageFilter.h"
#include "itkMultiplyImageFilter.h"

class landmarksConstellationDetector
{
  using MatrixType = vnl_matrix<double>;
  using VectorType = vnl_vector<double>;
  using ValMapType = std::map<std::string, float>;

public:
  landmarksConstellationDetector(const LandmarksMapType & orig_lmks)
    : m_mspQualityLevel(1)
    , m_orig_lmks_updated{ orig_lmks } // Initialize the updated landmarks
    , m_orig_lmks_forced{ orig_lmks }
    , m_HoughEyeFailure(false)
  {
    // Build midline landmarks name list
    this->m_MidlinePointsList.emplace_back("AC");
    this->m_MidlinePointsList.emplace_back("PC");
    this->m_MidlinePointsList.emplace_back("RP");
    this->m_MidlinePointsList.emplace_back("VN4");
    this->m_MidlinePointsList.emplace_back("aq_4V");
    this->m_MidlinePointsList.emplace_back("genu");
    this->m_MidlinePointsList.emplace_back("rostrum");
    this->m_MidlinePointsList.emplace_back("BPons");
    this->m_MidlinePointsList.emplace_back("optic_chiasm");
    this->m_MidlinePointsList.emplace_back("mid_ant");
    this->m_MidlinePointsList.emplace_back("mid_horiz");
    this->m_MidlinePointsList.emplace_back("mid_prim");
    this->m_MidlinePointsList.emplace_back("mid_prim_inf");
    this->m_MidlinePointsList.emplace_back("mid_prim_sup");
    this->m_MidlinePointsList.emplace_back("mid_sup");
  }

  const LandmarksMapType &
  Getorig_lmks_updated() const
  {
    return this->m_orig_lmks_updated;
  }

  void
  SetMSPQualityLevel(const int newLevel)
  {
    this->m_mspQualityLevel = newLevel;
  }

  void
  SetTemplateRadius(const ValMapType & NewTemplateRadius)
  {
    m_TemplateRadius = NewTemplateRadius;
  }

  void
  SeteyeFixed_img(SImageType::Pointer image)
  {
    m_eyeFixed_img = image; // This input is the output of Hough Eye detector
  }

  void
  SetInputTemplateModel(landmarksConstellationModelIO & myModel)
  {
    m_InputTemplateModel = myModel;
  }

  double
  GetModelRadius(std::string PointName)
  {
    return m_InputTemplateModel.GetRadius(PointName);
  }

  RigidTransformType::Pointer
  Getorig2msp_img_tfm() const
  {
    RigidTransformType::Pointer value = RigidTransformType::New();

    value->SetFixedParameters(this->m_eyeFixed2msp_img_tfm->GetFixedParameters());
    value->SetParameters(this->m_eyeFixed2msp_img_tfm->GetParameters());
    return value;
  }

  void
  SetResultsDir(std::string resultsDir)
  {
    this->m_ResultsDir = resultsDir;
  }

  void
  Setorig2eyeFixed_img_tfm(const VersorTransformType::Pointer houghEyeTransform)
  {
    this->m_orig2eyeFixed_img_tfm = houghEyeTransform;
  }

  void
  SetHoughEyeFailure(bool failure)
  {
    m_HoughEyeFailure = failure;
  }

  void
  SetLlsMeans(std::map<std::string, std::vector<double>> & llsMeans)
  {
    this->m_LlsMeans = llsMeans;
  }

  void
  SetLlsMatrices(std::map<std::string, MatrixType> & llsMatrices)
  {
    this->m_LlsMatrices = llsMatrices;
  }

  /** Set search radii for corresponding landmarks **/
  void
  SetSearchRadii(const std::map<std::string, double> & radii)
  {
    this->m_SearchRadii = radii;
  }

  void
  Compute(SImageType::Pointer orig_space_image);

  SImageType::Pointer
  GetTaggedImage(SImageType::Pointer original_space_image) const
  {
    itk::ImageDuplicator<SImageType>::Pointer duplicator = itk::ImageDuplicator<SImageType>::New();

    duplicator->SetInputImage(original_space_image);
    SImageType::Pointer taggedImage = duplicator->GetOutput();

    SImageType::PixelType low = 0;
    SImageType::PixelType high = 0;
    setLowHigh<SImageType>(taggedImage, low, high, 0.01F);

    SImageType::IndexType PTIndex;
    taggedImage->TransformPhysicalPointToIndex(GetNamedPointFromLandmarkList(this->Getorig_lmks_updated(), "AC"),
                                               PTIndex);
    taggedImage->SetPixel(PTIndex, high);
    taggedImage->TransformPhysicalPointToIndex(GetNamedPointFromLandmarkList(this->Getorig_lmks_updated(), "PC"),
                                               PTIndex);
    taggedImage->SetPixel(PTIndex, high);
    taggedImage->TransformPhysicalPointToIndex(GetNamedPointFromLandmarkList(this->Getorig_lmks_updated(), "VN4"),
                                               PTIndex);
    taggedImage->SetPixel(PTIndex, high);
    taggedImage->TransformPhysicalPointToIndex(GetNamedPointFromLandmarkList(this->Getorig_lmks_updated(), "RP"),
                                               PTIndex);
    taggedImage->SetPixel(PTIndex, high);
    return taggedImage;
  }

  void
  SetatlasVolume(const std::string & atlasVolume)
  {
    this->m_atlasVolume = atlasVolume;
  }
  void
  SetatlasLandmarks(const std::string & atlasLandmarks)
  {
    this->m_atlasLandmarks = atlasLandmarks;
  }
  void
  SetatlasLandmarkWeights(const std::string & atlasLandmarkWeights)
  {
    this->m_atlasLandmarkWeights = atlasLandmarkWeights;
  }

  VersorTransformType::Pointer
  GetImageOrigToACPCVersorTransform() const;
  void
  ComputeFinalRefinedACPCAlignedTransform(SImageType::Pointer      original_space_image,
                                          const LandmarksMapType & updated_orig_lmks);

protected:
private:
  void
  EulerToVersorRigid(VersorTransformType::Pointer &, const RigidTransformType::ConstPointer);

  void
  DoResampleInPlace(const SImageType::ConstPointer, const RigidTransformType::ConstPointer, SImageType::Pointer &);

  static VersorTransformType::Pointer
  Compute_orig2msp_img_tfm(const SImagePointType & RP, const SImagePointType & AC, const SImagePointType & PC);

  static bool
  mapHasKey(const LandmarksMapType & map, const std::string key)
  {
    return map.find(key) != map.cend();
  }

  static VersorTransformType::Pointer
  GetLandmarkTransformFromImageTransform(VersorTransformType::ConstPointer orig2msp_img_tfm);
  // Linear model estimation using EPCA
  void
  LinearEstimation(LandmarksMapType &               msp_lmks_linearly_estimated,
                   const std::vector<std::string> & processingList,
                   unsigned int                     numBasePoints);

  /*
   * Mainly designed for estimating the search center of
   * 4VN, ac, pc points in a morphometric way.
   *
   * Solve 2rd-order equations for BC.
   * Given:
   * BC.GetNorm() = BCMean.GetNorm();
   * angle(ABC) = angle(ABCMean);
   * BA[0] = BC[0] = 0.
   *
   * Then, we have
   * BA * BC = BA.GetNorm() * BC.GetNorm() * cosTheta;
   * Rearrange to solve a*BC[2]^2 + b*BC[2] + c = 0
   */
  SImageType::PointType::VectorType
  FindVectorFromPointAndVectors(SImageType::PointType::VectorType BA,
                                SImageType::PointType::VectorType BAMean,
                                SImageType::PointType::VectorType BCMean,
                                int                               sign);

  SImageType::PointType
  FindCandidatePoints(SImageType::Pointer volumeMSP,
                      SImageType::Pointer mask_LR,
                      const double        LR_restrictions,
                      const double        PA_restrictions,
                      const double        SI_restrictions,
                      // INFO: restrictions should really be ellipsoidal values
                      const SImageType::PointType::VectorType &                      CenterOfSearchArea,
                      const std::vector<std::vector<float>> &                        TemplateMean,
                      const landmarksConstellationModelIO::IndexLocationVectorType & model,
                      double &                                                       cc_Max,
                      const std::string &                                            mapID);

  SImageType::Pointer           m_eyeFixed_img;
  SImageType::Pointer           m_msp_img;
  landmarksConstellationModelIO m_InputTemplateModel;
  ValMapType                    m_TemplateRadius;
  int                           m_mspQualityLevel;
  std::string                   m_ResultsDir;


  LandmarksMapType       m_orig_lmks_updated;
  const LandmarksMapType m_orig_lmks_forced; // named points in the original space

  // name list of the landmarks that should be treated as midline landmarks
  std::vector<std::string> m_MidlinePointsList;

  RigidTransformType::Pointer  m_eyeFixed2msp_img_tfm;
  VersorTransformType::Pointer m_orig2msp_img_tfm;

  VersorTransformType::Pointer m_orig2eyeFixed_img_tfm;
  bool                         m_HoughEyeFailure;

  // Store linear model parameters
  // Note each matrix of m_LlsMatrices is actually cascaded by two mapping:
  // input space -> principal components space, and the coefs minimize their
  // least squares
  std::map<std::string, MatrixType> m_LlsMatrices;       // parameter
                                                         // matricess
  std::map<std::string, std::vector<double>> m_LlsMeans; // means to be
                                                         // subtracted by
                                                         // PCA
  std::map<std::string, double> m_SearchRadii;

  std::string m_atlasVolume;          // The reference atlas image
  std::string m_atlasLandmarks;       // The reference atlas landmarks
  std::string m_atlasLandmarkWeights; // The reference atlas landmark weights
};

void
WriteManualFixFiles(const std::string &      EMSP_Fiducial_file_name,
                    SImageType * const       mspVolume,
                    const std::string &      resultDir,
                    const LandmarksMapType & errorLmks,
                    const std::string &      failureMessage,
                    const bool               throwException);

#endif // __landmarksConstellationDetector__h
