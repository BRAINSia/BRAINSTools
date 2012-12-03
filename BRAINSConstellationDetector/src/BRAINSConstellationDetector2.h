/*
 * Author: Wei Lu
 * at Psychiatry Imaging Lab,
 * University of Iowa Health Care 2010
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * Redistributions of source code must retain the above copyright notice, this
 * list of conditions and the following disclaimer.
 *
 * Redistributions in binary form must reproduce the above copyright notice,
 * this list of conditions and the following disclaimer in the documentation
 * and/or other materials provided with the distribution.
 *
 * Neither the name of the Nathan Kline Institute nor the names of its
 * contributors may be used to endorse or promote products derived from this
 * software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 */

#ifndef __BRAINSConstellationDetector2_h
#define __BRAINSConstellationDetector2_h

#include "itkImageToImageFilter.h"

#include "landmarksConstellationDetector.h"
#include "itkTransformFileWriter.h"
#include "itkIntensityWindowingImageFilter.h"
#include "itkNearestNeighborInterpolateImageFunction.h"
#include "itkMultiplyImageFilter.h"
#include "landmarkIO.h"
#include "itkLargestForegroundFilledMaskImageFilter.h"
#include "ChopImageBelowLowerBound.h"
#include "itkResampleInPlaceImageFilter.h"
#include "itkImageDuplicator.h"

#include <string>
#include <vector>
#include <map>

#include "BRAINSMacro.h"

namespace itk
{
// Software Guide : BeginLatex
//
// \class BRAINSConstellationDetector2
// \brief ...
//
// This filter derives from the base class ImageToImageFilter
// Please refer to BRAINSConstellationDetector for details
//
// This filter takes two input filenames and a bunch of parameters:
// 1) input image filename
// 2) model filename
// 3) other optional parameters
//
// This filter outputs the following data
// 1) Aligned image with the detection of RP, AC, PC points
// 2) other optional parameters
//
// \ingroup ImageFeatureExtraction
// \todo Update the doxygen documentation!!!
//
// Software Guide : EndLatex

template <class TInputImage, class TOutputImage>
class ITK_EXPORT BRAINSConstellationDetector2 :
  public         ImageToImageFilter<SImageType, SImageType>
{
public:

  /** Standard ITK typedef */
  typedef BRAINSConstellationDetector2               Self;
  typedef ImageToImageFilter<SImageType, SImageType> Superclass;
  typedef SmartPointer<Self>                         Pointer;
  typedef SmartPointer<const Self>                   ConstPointer;

  itkStaticConstMacro(Dimension, unsigned int, SImageType::ImageDimension);
  typedef SImageType::PointType                              SImagePointType;
  typedef vnl_matrix<double>                                 MatrixType;
  typedef std::map<std::string, SImagePointType>             LandmarksMapType;
  typedef ImageDuplicator<SImageType>                        DuplicatorType;
  typedef ResampleInPlaceImageFilter<SImageType, SImageType> ResampleIPFilterType;
  typedef typename ResampleIPFilterType::Pointer             ResampleIPFilterPointer;

  /** Run-time type information (and related methods) */
  itkTypeMacro(BRAINSConstellationDetector2, ImageToImageFilter);

  /** Method for creation through the object factory */
  itkNewMacro(Self);

  /** Display */
  void PrintSelf(std::ostream & os, Indent indent) const;

  // Set Basic Inputs

  /** Set the filename of the output transform */
  itkSetMacro(Transform, std::string);

  /** Set the filename of the model file */
  itkSetMacro(InputTemplateModel, std::string);

  /** Set MSP quality level */
  itkSetMacro(MspQualityLevel, unsigned int);

  /** Set Otsu percentile threshold */
  itkSetMacro(OtsuPercentileThreshold, double);

  /** Set AC lower bound */
  itkSetMacro(AcLowerBound, double);

  /** Set Cut Out Head In Output volumes */
  itkSetMacro(CutOutHeadInOutputVolume, bool);

  /** Set rescale intensities */
  itkSetMacro(RescaleIntensities, bool);

  /** Set trim rescaled intensities */
  itkSetMacro(TrimRescaledIntensities, double);

  /** Set rescale intensities output range */
  VECTORitkSetMacro(RescaleIntensitiesOutputRange, std::vector<int> );

  /** Set the background fill value string */
  itkSetMacro(BackgroundFillValueString, std::string);

  /** Set use windowed sinc */
  itkSetMacro(InterpolationMode, std::string);

  /** Set Hough eye transform */
  itkSetObjectMacro(HoughEyeTransform, VersorTransformType);

  /** Set LE point */
  itkSetMacro(LEPoint, SImagePointType);

  /** Set RE point */
  itkSetMacro(REPoint, SImagePointType);

  /** Set center of head mass **/
  itkSetMacro(CenterOfHeadMass, SImagePointType);

  /** Set the original input image before the Hough eye detector */
  itkSetObjectMacro(OriginalInputImage, SImageType);
  itkGetConstObjectMacro(OriginalInputImage, SImageType);

  // Get Basic Outputs

  /** Get the versor transform */
  itkGetObjectMacro(VersorTransform, VersorTransformType);
  itkGetConstObjectMacro(VersorTransform, VersorTransformType);
  itkGetObjectMacro(InvVersorTransform, VersorTransformType);
  itkGetConstObjectMacro(InvVersorTransform, VersorTransformType);

  /** Get the named points in original space */
  const LandmarksMapType & GetOriginalPoints()
  {
    return this->m_OriginalPoints;
  }

  /** Get the aligned named points */
  const LandmarksMapType & GetAlignedPoints()
  {
    return this->m_AlignedPoints;
  }

  // This is wrong?
  const SImagePointType & GetAlignedPoint(std::string name)
  {
    return this->m_AlignedPoints[name];
  }

  /** Get the interpolated output isotropic image */
  itkGetConstObjectMacro(OutputResampledImage, SImageType);

  /** Get the output untransformed clipped volume */
  itkGetConstObjectMacro(OutputUntransformedClippedVolume, SImageType);

  /** Get the Hough eye transform */
  itkGetObjectMacro(HoughEyeTransform, VersorTransformType);
  itkGetConstObjectMacro(HoughEyeTransform, VersorTransformType);

  /** Set the Hough eye failure report */
  void SetHoughEyeFailure(const bool failure)
  {
    this->m_HoughEyeFailure = failure;
  }

  /** Set llsMeans **/
  void SetLlsMeans(const std::map<std::string, std::vector<double> > & llsMeans)
  {
    this->m_LlsMeans = llsMeans;
  }

  /** Set llsMatrices **/
  void SetLlsMatrices(const std::map<std::string, MatrixType> & llsMatrices)
  {
    this->m_LlsMatrices = llsMatrices;
  }

  /** Set search radii for corresponding landmarks **/
  void SetSearchRadii(const std::map<std::string, double> & radii)
  {
    this->m_SearchRadii = radii;
  }

  /** Set AC mean **/
  itkSetMacro(ACMean, SImagePointType);

  // Set Advanced Inputs

  /** Set force AC point */
  VECTORitkSetMacro(ForceACPoint, std::vector<float> );

  /** Set force PC point */
  VECTORitkSetMacro(ForcePCPoint, std::vector<float> );

  /** Set force VN4 point */
  VECTORitkSetMacro(ForceVN4Point, std::vector<float> );

  /** Set force MPJ point */
  VECTORitkSetMacro(ForceRPPoint, std::vector<float> );

  /** Set MPJ search radius */
  itkSetMacro(RadiusMPJ, double);

  /** Set AC search radius */
  itkSetMacro(RadiusAC, double);

  /** Set PC search radius */
  itkSetMacro(RadiusPC, double);

  /** Set VN4 search radius */
  itkSetMacro(RadiusVN4, double);

  /** Set use debug mode */
  itkSetMacro(Debug, bool);

  /** Set use verbose mode */
  itkSetMacro(Verbose, bool);

  /** Set debugging images level */
  itkSetMacro(WritedebuggingImagesLevel, unsigned int);

  /** Set branded 2D image filename */
  itkSetMacro(WriteBranded2DImage, std::string);

  /** Set results dir */
  itkSetMacro(ResultsDir, std::string);

  /** Set/Get EMSP landmarks */
  void SetLandmarksEMSP(LandmarksMapType landmarks)
  {
    m_landmarksEMSP.clear();
    m_landmarksEMSP.insert( landmarks.begin(), landmarks.end() );
  }

  LandmarksMapType GetLandmarksEMSP()
  {
    return m_landmarksEMSP;
  }

protected:

  BRAINSConstellationDetector2();

  void GenerateData();

  /** Essential Parameters */
  // Inputs
  std::string      m_Transform;
  std::string      m_InputTemplateModel;
  unsigned int     m_MspQualityLevel;               // default = 2
  double           m_OtsuPercentileThreshold;       // default = 0.01
  double           m_AcLowerBound;                  // default = 1000.0
  bool             m_CutOutHeadInOutputVolume;      // default = false
  bool             m_RescaleIntensities;            // default = false
  double           m_TrimRescaledIntensities;       // default = 4.4172
  std::vector<int> m_RescaleIntensitiesOutputRange; // default = [40,
                                                    // 4000]
  std::string m_BackgroundFillValueString;          // default = "0"
  std::string m_InterpolationMode;                  // default = "Linear"

  // a local editable copy of original input before Hough eye detector
  // Note: this->GetInput() will return a const input after Hough eye.
  SImageType::Pointer m_OriginalInputImage;

  VersorTransformType::Pointer m_HoughEyeTransform; // help to get the points
                                                    // location in the original
                                                    // space
  SImagePointType m_LEPoint;                        // automated estimated LE in
                                                    // the original space by
                                                    // Hough eye detector
  SImagePointType  m_REPoint;
  SImagePointType  m_CenterOfHeadMass;
  LandmarksMapType m_landmarksEMSP;
  bool             m_HoughEyeFailure;
  std::map<std::string,
           MatrixType> m_LlsMatrices;
  std::map<std::string,
           std::vector<double> > m_LlsMeans;
  SImagePointType                m_ACMean;
  std::map<std::string, double>  m_SearchRadii;

  // Outputs
  VersorTransformType::Pointer m_VersorTransform;
  VersorTransformType::Pointer m_InvVersorTransform;
  LandmarksMapType             m_AlignedPoints;
  LandmarksMapType             m_OriginalPoints;
  SImageType::Pointer          m_OutputImage;          // Output image w/o
                                                       // interpolation
  SImageType::Pointer m_OutputResampledImage;          // Output image w/
                                                       // interpolation
  SImageType::Pointer m_OutputUntransformedClippedVolume;

  /** Advanced parameters */
  /** Manual Override */
  // Inputs
  std::vector<float> m_ForceACPoint;              // default = 0.
  std::vector<float> m_ForcePCPoint;              // default = 0.
  std::vector<float> m_ForceVN4Point;             // default = 0.
  std::vector<float> m_ForceRPPoint;              // default = 0.

  /** Model Override */
  // Inputs
  double m_RadiusMPJ;                         // default = -1.
  double m_RadiusAC;                          // default = -1.
  double m_RadiusPC;                          // default = -1.
  double m_RadiusVN4;                         // default = -1.

  /** Debug Options */
  // Inputs
  bool         m_Debug;                     // default = false
  bool         m_Verbose;                   // default = false
  unsigned int m_WritedebuggingImagesLevel; // default = 0
  std::string  m_WriteBranded2DImage;
  std::string  m_ResultsDir;                    // default = "./"
private:

  BRAINSConstellationDetector2(const Self &)
  {
  }

  void operator=(const Self &)
  {
  }
};
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "BRAINSConstellationDetector2.hxx"
#endif

#endif
