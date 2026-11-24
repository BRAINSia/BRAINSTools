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
 * Portions of the keypoint detector and descriptor are adapted from SIFT3D,
 * https://github.com/bbrister/SIFT3D, under the MIT License:
 *
 *   Copyright (c) 2015-2016 Blaine Rister et al.
 *
 *   Permission is hereby granted, free of charge, to any person obtaining a
 *   copy of this software and associated documentation files (the
 *   "Software"), to deal in the Software without restriction, including
 *   without limitation the rights to use, copy, modify, merge, publish,
 *   distribute, sublicense, and/or sell copies of the Software, and to permit
 *   persons to whom the Software is furnished to do so, subject to the
 *   following conditions:
 *
 *   The above copyright notice and this permission notice shall be included
 *   in all copies or substantial portions of the Software.
 *
 *   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
 *   OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 *   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN
 *   NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
 *   DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
 *   OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE
 *   USE OR OTHER DEALINGS IN THE SOFTWARE.
 *
 * Reference: B. Rister, M. A. Horowitz, and D. L. Rubin, "Volumetric Image
 * Registration From Invariant Keypoints," IEEE Trans. Image Process. 26(10),
 * 2017.
 */
#ifndef BRAINSSIFT3DFeatureDetector_h
#define BRAINSSIFT3DFeatureDetector_h

#include "itkImage.h"
#include "itkObject.h"
#include "itkObjectFactory.h"
#include "itkPoint.h"

#include <vector>

namespace brains
{

/** Keypoints and 768-element descriptors of one image, in physical (LPS) space. */
struct SIFT3DFeatures
{
  static constexpr unsigned int DescriptorLength = 768;
  using PointType = itk::Point<double, 3>;

  std::vector<PointType> locations;
  std::vector<double>    scales;      // mm
  std::vector<double>    strengths;   // |DoG| at the extremum
  std::vector<float>     descriptors; // row-major, size() x DescriptorLength

  size_t
  size() const
  {
    return locations.size();
  }
};

/** \class SIFT3DFeatureDetector
 *
 * 3D scale-invariant keypoints with rotation-invariant descriptors, after
 * SIFT3D (Rister et al. 2017): difference-of-Gaussian extrema, orientation
 * from the eigenvectors of the windowed structure tensor with a corner test,
 * and a 4x4x4 grid of 12-bin icosahedral gradient histograms.
 *
 * The input is first resampled onto an isotropic grid with identity
 * direction, so voxel order is the same for every input (descriptors are not
 * mirror-invariant) and all lengths are in grid voxels. The scale space is
 * the standard one: each octave is blurred in its own voxel units and
 * edges are clamped.
 */
class SIFT3DFeatureDetector : public itk::Object
{
public:
  ITK_DISALLOW_COPY_AND_MOVE(SIFT3DFeatureDetector);

  using Self = SIFT3DFeatureDetector;
  using Superclass = itk::Object;
  using Pointer = itk::SmartPointer<Self>;
  using ConstPointer = itk::SmartPointer<const Self>;
  using ImageType = itk::Image<float, 3>;

  itkNewMacro(Self);
  itkOverrideGetNameOfClassMacro(SIFT3DFeatureDetector);

  /** Isotropic spacing (mm) of the grid features are computed on. */
  itkSetMacro(SamplingSpacing, double);
  itkGetConstMacro(SamplingSpacing, double);
  /** DoG extrema below this fraction of the level's largest |DoG| are dropped. */
  itkSetMacro(PeakThreshold, double);
  itkGetConstMacro(PeakThreshold, double);
  /** Minimum corner score (|cos| between mean gradient and eigenvectors). */
  itkSetMacro(CornerThreshold, double);
  itkGetConstMacro(CornerThreshold, double);
  itkSetMacro(NumberOfKeypointLevels, unsigned int);
  itkGetConstMacro(NumberOfKeypointLevels, unsigned int);
  itkSetMacro(Sigma0, double);
  itkGetConstMacro(Sigma0, double);
  itkSetMacro(SigmaN, double);
  itkGetConstMacro(SigmaN, double);
  /** Invert contrast inside the head (p99 - I) and zero the background, so a
   *  T1-weighted image matches T2-weighted edge polarity. */
  itkSetMacro(InvertIntensity, bool);
  itkGetConstMacro(InvertIntensity, bool);
  itkBooleanMacro(InvertIntensity);
  /** Keep only the strongest keypoints; 0 keeps all. */
  itkSetMacro(MaximumNumberOfKeypoints, unsigned int);
  itkGetConstMacro(MaximumNumberOfKeypoints, unsigned int);

  SIFT3DFeatures
  Compute(const ImageType * image) const;

protected:
  SIFT3DFeatureDetector() = default;
  ~SIFT3DFeatureDetector() override = default;
  void
  PrintSelf(std::ostream & os, itk::Indent indent) const override;

private:
  double       m_SamplingSpacing{ 1.0 };
  double       m_PeakThreshold{ 0.1 };
  double       m_CornerThreshold{ 0.4 };
  unsigned int m_NumberOfKeypointLevels{ 3 };
  double       m_Sigma0{ 1.6 };
  double       m_SigmaN{ 1.15 };
  unsigned int m_MaximumNumberOfKeypoints{ 0 };
  bool         m_InvertIntensity{ false };
};

} // namespace brains

#endif
