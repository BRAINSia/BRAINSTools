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
 * Author: Han J. Johnson, Wei Lu
 * at Psychiatry Imaging Lab,
 * University of Iowa Health Care 2010
 */

#ifndef __landmarksConstellationDetector__h
#define __landmarksConstellationDetector__h

#include "landmarksConstellationCommon.h"
#include "landmarksConstellationModelIO.h"
#include "itkAffineTransform.h"
#include "itkOtsuThresholdImageFilter.h"
#include "Slicer3LandmarkIO.h"
#include "PrepareOutputImages.h"

#include <map>

#include "itkMaskedFFTNormalizedCorrelationImageFilter.h"
#include "itkMinimumMaximumImageCalculator.h"
#include "itkBinaryImageToLabelMapFilter.h"
#include "itkLabelMapToLabelImageFilter.h"
#include "itkLabelStatisticsImageFilter.h"
#include "itkMultiplyImageFilter.h"

class landmarksConstellationDetector
{
  typedef vnl_matrix<double>                           MatrixType;
  typedef vnl_vector<double>                           VectorType;
  typedef std::map<std::string, float>                 ValMapType;
public:
  landmarksConstellationDetector() :
    m_mspQualityLevel(1),
    m_HoughEyeFailure(false)
  {
    // Build midline landmarks name list
    this->m_MidlinePointsList.push_back("AC");
    this->m_MidlinePointsList.push_back("PC");
    this->m_MidlinePointsList.push_back("RP");
    this->m_MidlinePointsList.push_back("VN4");
    this->m_MidlinePointsList.push_back("aq_4V");
    this->m_MidlinePointsList.push_back("genu");
    this->m_MidlinePointsList.push_back("rostrum");
    this->m_MidlinePointsList.push_back("BPons");
    this->m_MidlinePointsList.push_back("optic_chiasm");
    this->m_MidlinePointsList.push_back("mid_ant");
    this->m_MidlinePointsList.push_back("mid_horiz");
    this->m_MidlinePointsList.push_back("mid_prim");
    this->m_MidlinePointsList.push_back("mid_prim_inf");
    this->m_MidlinePointsList.push_back("mid_prim_sup");
    this->m_MidlinePointsList.push_back("mid_sup");
  }

  // Force the setting of the point values to override those that were specified.
  void SetOriginalSpaceNamedPoint(const std::string & NamedPoint, const SImageType::PointType & PointValue)
  {
    this->m_OriginalSpaceNamedPoints[NamedPoint] = PointValue;
    return;
  }

  const LandmarksMapType & GetOriginalSpaceNamedPoints(void) const
  {
    return this->m_OriginalSpaceNamedPoints;
  }

  void SetMSPQualityLevel(const int newLevel)
  {
    this->m_mspQualityLevel = newLevel;
  }

  void SetTemplateRadius(const ValMapType & NewTemplateRadius)
  {
    m_TemplateRadius = NewTemplateRadius;
  }

  void SetVolumeRoughAlignedWithHoughEye(SImageType::Pointer image)
  {
    m_VolumeRoughAlignedWithHoughEye = image;  // This input is the output of Hough Eye detector
  }

  void SetOriginalInputImage(SImageType::Pointer image)
  {
    m_OriginalInputImage = image; // This is the original input of BCD
  }

  SImageType::Pointer GetOriginalInputImage(void) const
  {
    return this->m_OriginalInputImage; // Returns the original input of BCD
  }

  void SetInputTemplateModel(landmarksConstellationModelIO & myModel)
  {
    m_InputTemplateModel = myModel;
  }

  double GetModelHeight(std::string PointName)
  {
    return m_InputTemplateModel.GetHeight(PointName);
  }

  double GetModelRadius(std::string PointName)
  {
    return m_InputTemplateModel.GetRadius(PointName);
  }

  RigidTransformType::Pointer GetTransformToMSP(void) const
  {
    RigidTransformType::Pointer value = RigidTransformType::New();

    value->SetFixedParameters( this->m_finalTmsp->GetFixedParameters() );
    value->SetParameters( this->m_finalTmsp->GetParameters() );
    return value;
  }

  SImageType::PointType GetCenterOfHeadMassMSP(void) const
  {
    return this->m_CenterOfHeadMassEMSP;
  }

  void SetResultsDir(std::string resultsDir)
  {
    this->m_ResultsDir = resultsDir;
  }

  void SetHoughEyeTransform(VersorTransformType::Pointer houghEyeTransform)
  {
    this->m_HoughEyeTransform = houghEyeTransform;
  }

  void SetHoughEyeFailure(bool failure)
  {
    m_HoughEyeFailure = failure;
  }

  void SetLEPoint(const SImageType::PointType & LEPoint)
  {
    this->m_LEPoint = LEPoint;
  }

  const SImageType::PointType & GetLEPoint() const
  {
    return this->m_LEPoint;
  }

  void SetREPoint(const SImageType::PointType & REPoint)
  {
    this->m_REPoint = REPoint;
  }

  const SImageType::PointType & GetREPoint() const
  {
    return this->m_REPoint;
  }

  void SetCenterOfHeadMass(const SImageType::PointType & centerOfHeadMass)
  {
    m_CenterOfHeadMass = centerOfHeadMass;
  }

  void SetLlsMeans(std::map<std::string, std::vector<double> > & llsMeans)
  {
    this->m_LlsMeans = llsMeans;
  }

  void SetLlsMatrices(std::map<std::string, MatrixType> & llsMatrices)
  {
    this->m_LlsMatrices = llsMatrices;
  }

  /** Set search radii for corresponding landmarks **/
  void SetSearchRadii(const std::map<std::string, double> & radii)
  {
    this->m_SearchRadii = radii;
  }

  void SetLandmarksEMSP(LandmarksMapType landmarks)
  {
    m_NamedPointEMSP.clear();
    m_NamedPointEMSP.insert( landmarks.begin(), landmarks.end() );
  }

  void Compute(void);

  SImageType::Pointer GetTaggedImage(void) const
  {
    itk::ImageDuplicator<SImageType>::Pointer duplicator = itk::ImageDuplicator<SImageType>::New();

    duplicator->SetInputImage(this->GetOriginalInputImage());
    SImageType::Pointer taggedImage = duplicator->GetModifiableOutput();

    SImageType::PixelType low=0.0;
    SImageType::PixelType high=0.0;
    setLowHigh<SImageType>(taggedImage, low, high, 0.01F);

    SImageType::IndexType PTIndex;
    taggedImage->TransformPhysicalPointToIndex(GetNamedPointFromLandmarkList(this->GetOriginalSpaceNamedPoints(),"AC"), PTIndex);
    taggedImage->SetPixel(PTIndex, high);
    taggedImage->TransformPhysicalPointToIndex(GetNamedPointFromLandmarkList(this->GetOriginalSpaceNamedPoints(),"PC"), PTIndex);
    taggedImage->SetPixel(PTIndex, high);
    taggedImage->TransformPhysicalPointToIndex(GetNamedPointFromLandmarkList(this->GetOriginalSpaceNamedPoints(),"VN4"), PTIndex);
    taggedImage->SetPixel(PTIndex, high);
    taggedImage->TransformPhysicalPointToIndex(GetNamedPointFromLandmarkList(this->GetOriginalSpaceNamedPoints(),"RP"), PTIndex);
    taggedImage->SetPixel(PTIndex, high);
    return taggedImage;
  }


  SImageType::Pointer GetVolumeMSP()
  {
    return this->m_VolumeMSP;
  }

  void SetatlasVolume( const std::string & atlasVolume )
    {
    this->m_atlasVolume = atlasVolume;
    }
  void SetatlasLandmarks( const std::string & atlasLandmarks )
    {
    this->m_atlasLandmarks = atlasLandmarks;
    }
  void SetatlasLandmarkWeights( const std::string & atlasLandmarkWeights )
    {
    this->m_atlasLandmarkWeights = atlasLandmarkWeights;
    }

  VersorTransformType::Pointer GetImageOrigToACPCVersorTransform(void) const;
  void ComputeFinalRefinedACPCAlignedTransform(void);
protected:

private:
  void EulerToVersorRigid( VersorTransformType::Pointer &, const RigidTransformType::ConstPointer );

  void DoResampleInPlace( const SImageType::ConstPointer, const RigidTransformType::ConstPointer, SImageType::Pointer & );

  VersorTransformType::Pointer ComputeACPCAlignedZeroCenteredTransform(void);

  // Linear model estimation using EPCA
  void LinearEstimation( LandmarksMapType & namedPoints, const std::vector<std::string> & processingList,
                         unsigned int numBasePoints );

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
  SImageType::PointType::VectorType FindVectorFromPointAndVectors(SImageType::PointType::VectorType BA,
                                                                  SImageType::PointType::VectorType BAMean,
                                                                  SImageType::PointType::VectorType BCMean, int sign);

  SImageType::PointType FindCandidatePoints(SImageType::Pointer volumeMSP, SImageType::Pointer mask_LR,
                                            const double LR_restrictions, const double PA_restrictions,
                                            const double SI_restrictions,
                                            // TODO: restrictions should really be ellipsoidal values
                                            const SImageType::PointType::VectorType & CenterOfSearchArea,
                                            const std::vector<std::vector<float> > & TemplateMean,
                                            const landmarksConstellationModelIO::IndexLocationVectorType & model,
                                            double & cc_Max,
                                            const std::string & mapID);

  RigidTransformType::Pointer   m_TmspBasedOnReflectionCrossCorrelation;
  SImageType::PointType         m_CenterOfHeadMassEMSP;
  SImageType::PointType         m_BestCenter;
  SImageType::Pointer           m_VolumeRoughAlignedWithHoughEye;
  SImageType::Pointer           m_OriginalInputImage;
  SImageType::Pointer           m_VolumeMSP;
  landmarksConstellationModelIO m_InputTemplateModel;
  ValMapType                    m_TemplateRadius;
  int                           m_mspQualityLevel;
  std::string                   m_ResultsDir;

  LandmarksMapType m_OriginalSpaceNamedPoints;             // named points in the
                                                  // original space
                                                  // even before the Hough eye
                                                  // detector

  LandmarksMapType m_NamedPointEMSP;              // named points in EMSP space

  std::vector<std::string> m_MidlinePointsList;   // name list of the landmarks
                                                  // that
                                                  // should be treated as
                                                  // midline landmarks

  RigidTransformType::Pointer  m_finalTmsp;
  VersorTransformType::Pointer m_ImageOrigToACPCVersorTransform;

  // Wei: Read in LE, RE value for linear model estimation
  VersorTransformType::Pointer m_HoughEyeTransform;
  bool                         m_HoughEyeFailure;
  SImageType::PointType        m_LEPoint;    // in input space
  SImageType::PointType        m_REPoint;

  SImageType::PointType m_ReferencePointAC;
  SImageType::PointType m_ReferencePointPC;
  SImageType::PointType m_CenterOfHeadMass;

  // Store linear model parameters
  // Note each matrix of m_LlsMatrices is actually cascaded by two mapping:
  // input space -> principal components space, and the coefs minimize their
  // least squares
  std::map<std::string, MatrixType> m_LlsMatrices;              // parameter
                                                                // matricess
  std::map<std::string, std::vector<double> > m_LlsMeans;       // means to be
                                                                // subtracted by
                                                                // PCA
  std::map<std::string, double> m_SearchRadii;

  std::string m_atlasVolume; // The reference atlas image
  std::string m_atlasLandmarks; // The reference atlas landmarks
  std::string m_atlasLandmarkWeights; // The reference atlas landmark weights

};

//TODO:  Move out of class all together
void WriteManualFixFiles(const std::string &EMSP_Fiducial_file_name, SImageType * const mspVolume,
                         const std::string &resultDir, const LandmarksMapType & errorLmks,
                         const std::string & failureMessage) ;

#endif // __landmarksConstellationDetector__h
