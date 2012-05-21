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

#include <map>

class landmarksConstellationDetector
{
  typedef vnl_matrix<double>                           MatrixType;
  typedef vnl_vector<double>                           VectorType;
  typedef std::map<std::string, SImageType::PointType> LandmarksMapType;
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

  /**
   * Returns the named point (AC,PC,VN4,RP) in the space of the
   * original image.
   *
   * @author hjohnson (12/14/2008)
   *
   * @param NamedPoint
   *
   * @return SImageType::PointType
   */
  SImageType::PointType GetOriginalSpaceNamedPoint(const std::string & NamedPoint) const
  {
    LandmarksMapType::const_iterator itpair = this->m_NamedPoint.find(NamedPoint);

    if( itpair == this->m_NamedPoint.end() )
      {
      std::cout << "ERROR:  " << NamedPoint << " not found in list." << std::endl;
      return SImageType::PointType();
      }
    return itpair->second;
  }

  // Force the setting of the point values to override those that were
  // specified.
  void SetOriginalSpaceNamedPoint(const std::string & NamedPoint, const SImageType::PointType & PointValue)
  {
    this->m_NamedPoint[NamedPoint] = PointValue;
    return;
  }

  const LandmarksMapType & GetNamedPoints()
  {
    return this->m_NamedPoint;
  }

  void SetMSPQualityLevel(const int newLevel)
  {
    this->m_mspQualityLevel = newLevel;
  }

  void SetTemplateRadius(const ValMapType & NewTemplateRadius)
  {
    m_TemplateRadius = NewTemplateRadius;
  }

  void SetVolOrig(SImageType::Pointer image)
  {
    m_VolOrig = image;  // This input is the output of Hough Eye detector
  }

  void SetOriginalInput(SImageType::Pointer image)
  {
    m_OriginalInput = image; // This is the original input of BCD
  }

  SImageType::Pointer GetOriginalInput(void) const
  {
    return this->m_OriginalInput; // Returns the original input of BCD
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

    duplicator->SetInputImage(this->m_VolOrig);
    duplicator->Update();
    SImageType::Pointer taggedImage = duplicator->GetOutput();

    SImageType::PixelType low, high;
    setLowHigh<SImageType>(taggedImage, low, high, 0.01F);

    SImageType::IndexType PTIndex;
    taggedImage->TransformPhysicalPointToIndex(this->GetOriginalSpaceNamedPoint("AC"), PTIndex);
    taggedImage->SetPixel(PTIndex, high);
    taggedImage->TransformPhysicalPointToIndex(this->GetOriginalSpaceNamedPoint("PC"), PTIndex);
    taggedImage->SetPixel(PTIndex, high);
    taggedImage->TransformPhysicalPointToIndex(this->GetOriginalSpaceNamedPoint("VN4"), PTIndex);
    taggedImage->SetPixel(PTIndex, high);
    taggedImage->TransformPhysicalPointToIndex(this->GetOriginalSpaceNamedPoint("RP"), PTIndex);
    taggedImage->SetPixel(PTIndex, high);
    return taggedImage;
  }

  RigidTransformType::Pointer GetACPCAlignedZeroCenteredTransform(void) const
  {
    SImageType::PointType ZeroCenter;

    ZeroCenter.Fill(0.0);
    RigidTransformType::Pointer
      landmarkDefinedACPCAlignedToZeroTransform =
      computeTmspFromPoints(this->GetOriginalSpaceNamedPoint("RP"),
                            this->GetOriginalSpaceNamedPoint("AC"),
                            this->GetOriginalSpaceNamedPoint("PC"),
                            ZeroCenter);
    return landmarkDefinedACPCAlignedToZeroTransform;
  }

  SImageType::Pointer GetVolumeMSP()
  {
    return this->m_VolumeMSP;
  }

private:

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
                                            const bool ComputeOutsideSearchRadius, double & cc_Max,
                                            const std::string & mapID);

private:
  RigidTransformType::Pointer   m_TmspBasedOnReflectionCrossCorrelation;
  SImageType::PointType         m_CenterOfHeadMassEMSP;
  SImageType::PointType         m_BestCenter;
  SImageType::Pointer           m_VolOrig;
  SImageType::Pointer           m_OriginalInput;
  SImageType::Pointer           m_VolumeMSP;
  landmarksConstellationModelIO m_InputTemplateModel;
  ValMapType                    m_TemplateRadius;
  int                           m_mspQualityLevel;
  std::string                   m_ResultsDir;

  LandmarksMapType m_NamedPoint;                  // named points in the
                                                  // original space
                                                  // even before the Hough eye
                                                  // detector
  LandmarksMapType m_NamedPointEMSP;              // named points in EMSP space
  LandmarksMapType m_NamedPointACPC;              // named points in
                                                  // ACPC-aligned space
  LandmarksMapType m_NamedPointACPCRaw;           // reserve this map to avoid
                                                  // the
                                                  // accumulation of local
                                                  // search
                                                  // errors
  std::vector<std::string> m_MidlinePointsList;   // name list of the landmarks
                                                  // that
                                                  // should be treated as
                                                  // midline landmarks

  RigidTransformType::Pointer  m_finalTmsp;
  VersorTransformType::Pointer m_VersorTransform;
  VersorTransformType::Pointer m_InvVersorTransform;

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
};

#endif // __landmarksConstellationDetector__h
