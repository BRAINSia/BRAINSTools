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

/**
 * This is a temporary file to be used during refactoring to help consolidate a family of functions into a common
 * library.
 */
#ifndef __landmarksConstellationCommon_h
#define __landmarksConstellationCommon_h

#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <algorithm>

#include <numeric>
#include <itksys/SystemTools.hxx>
#include <string>

#include "itkImage.h"
#include "itkIO.h"

#include <itkThresholdImageFilter.h>
#include <itkImageDuplicator.h>
#include "itkMultiResolutionPyramidImageFilter.h"
#include <itkDiscreteGaussianImageFilter.h>
#include <itkMersenneTwisterRandomVariateGenerator.h>
#include <itkImageFileWriter.h>
#include <itkSimilarity3DTransform.h>
#include <itkVersorRigid3DTransform.h>

#include "itkSpatialObjectToImageFilter.h"
#include "itkRGBPixel.h"
#include "itkExtractImageFilter.h"
#include "itkStatisticsImageFilter.h"
#include "itkFlipImageFilter.h"

#include "itkImageDuplicator.h"
#include "itkEuler3DTransform.h"
#include <vnl/vnl_cross.h>

#include <itkScalarImageToHistogramGenerator.h>

#include "GenericTransformImage.h"
#include "itkLandmarkBasedTransformInitializer.h"
#include "itkSimilarity3DTransform.h"

#include "Slicer3LandmarkIO.h"

using ImagePointType = itk::Point<double, 3>;
using PointList = std::vector<ImagePointType>;
using VersorRigidTransformType = itk::VersorRigid3DTransform<double>;
using SimilarityTransformType = itk::Similarity3DTransform<double>;
using Euler3DTransformType = itk::Euler3DTransform<double>;

extern const unsigned int NO;
namespace LMC
{
extern bool debug;
extern bool globalverboseFlag;
} // namespace LMC

//
//
// ////////////////////////////////////////////////////////////////////////////////////////////////

namespace // avoid 'shadows declaration' warnings.
{

using SImageType = itk::Image<short, 3>;
using DImageType3D = itk::Image<double, 3>;
using FImageType3D = itk::Image<float, 3>;

using ByteImageType = itk::Image<unsigned char, 3>;

using SImagePointType = SImageType::PointType;

using RGBPixelType = itk::RGBPixel<unsigned char>;
using RGBImageType = itk::Image<RGBPixelType, 3>;
using RGB2DImageType = itk::Image<RGBPixelType, 2>;



using PyramidFilterType = itk::MultiResolutionPyramidImageFilter<SImageType, SImageType>;
using LinearInterpolatorType = itk::LinearInterpolateImageFunction<SImageType, double>;
} // namespace

#include "landmarksConstellationModelIO.h"

extern std::string globalResultsDir;
extern int         globalImagedebugLevel;
extern PyramidFilterType::Pointer
MakeThreeLevelPyramid(const SImageType::Pointer & refImage);

extern PyramidFilterType::Pointer
MakeOneLevelPyramid(const SImageType::Pointer & refImage);

extern SImageType::PointType
GetImageCenterPhysicalPoint(SImageType::Pointer & image);

extern unsigned char
ShortToUChar(short in, short min, short max);

extern SImageType::Pointer
CreateTestCenteredRotatedImage2(const VersorRigidTransformType::Pointer & ACPC_MSP_AlignedTransform,
                                /* const
                                  SImageType::PointType
                                  finalPoint, */
                                const SImageType::PointType         PreMSP_Point,
                                /*const*/ SImageType::Pointer &     image,
                                const Euler3DTransformType::Pointer & Point_Rotate);

extern void
decomposeRPAC(const SImageType::PointType & RP,
              const SImageType::PointType & PC,
              const SImageType::PointType & AC,
              double * const                RPPC_to_RPAC_angle,
              double * const                RPAC_over_RPPC);

extern void
MakeLabelImage(const SImageType::Pointer &   in,
               const SImageType::PointType & RP,
               const SImageType::PointType & AC,
               const SImageType::PointType & PC,
               const SImageType::PointType & VN4,
               const std::string &           fname);

extern VersorRigidTransformType::Pointer
ComputeMSP(SImageType::Pointer           input_image,
           const SImageType::PointType & input_image_center_of_mass,
           const int                     qualityLevel,
           double &                      cc);

extern SImageType::Pointer
CreatedebugPlaneImage(SImageType::Pointer                 referenceImage,
                      const VersorRigidTransformType::Pointer & MSPTransform,
                      const std::string &                 debugfilename);

extern void
CreatedebugPlaneImage(SImageType::Pointer referenceImage, const std::string & debugfilename);

extern VersorRigidTransformType::Pointer
computeTmspFromPoints_Versor(SImageType::PointType RP,
                      SImageType::PointType AC,
                      SImageType::PointType PC,
                      SImageType::PointType DesiredCenter);

extern Euler3DTransformType::Pointer
computeTmspFromPoints(SImageType::PointType RP,
                      SImageType::PointType AC,
                      SImageType::PointType PC,
                      SImageType::PointType DesiredCenter);

extern VersorRigidTransformType::Pointer
GetACPCAlignedZeroCenteredTransform(const LandmarksMapType & landmarks);

extern SImageType::PointType
GetCenterOfHeadMass(SImageType::Pointer volume);

extern SImageType::Pointer
MakeIsoTropicReferenceImage();

/*****************************************************************************/
template <typename ValuesType>
ValuesType
vectorNorm(const std::vector<ValuesType> & x)
{
  ValuesType norm = 0.0;

  for (auto it = x.begin(); it != x.end(); ++it)
  {
    const ValuesType value = *it;
    norm += value * value;
  }
  return std::sqrt(norm);
}

/*
 This function takes an input array (y) of size n and computes its L2
 norm.  Then, the array y is modified by dividing each of its elements by
 the L2 norm.  The resulting array will have an L2 norm of 1.0.

 Example:
 {
 float y[3]={0.0, 4.0, 3.0}
 double mean;

 normalizeVector(y,3);

 printf("y = {%f, %f, %f}\n", y[0],y[1],y[2]);
 }

 In this example, the L2 norm of y is 5.0.  Array y is replaced by: {0.0,
 0.8, 0.6}, which has an L2 norm of 1.0.

 If we run removeVectorMean(y,3) before calling normalizeVector(y,3), we
 will obtain an array which has a zero mean as well as unit norm.
 */
template <typename ValuesType>
void
normalizeVector(std::vector<ValuesType> & x)
{
  const ValuesType norm = vectorNorm(x);

  if (norm < std::numeric_limits<ValuesType>::epsilon())
  {
    std::cout << "WARNING:  ZERO NORM VECTOR." << __FILE__ << __LINE__ << std::endl;
    return;
  }
  for (auto it = x.begin(); it != x.end(); ++it)
  {
    *it /= norm;
  }
}

// ///////////////////////////////////////////////////////
// y is a (nxp) matrix.  This function operates on the columns of y.
// Upon completion, each of the p columns of y will have zero mean.
/*
 This function takes an input array (y) of size n, computes and returns
 the mean of the array.  Also, the input array y is modified by
 subtracting the computed mean from each of its elements.

 Example:

 {
 float y[3]={1.0, 5.0, 3.0}
 double mean;

 mean = removeVectorMean(y,3);

 printf("mean = %lf\n", mean);
 printf("y = {%f, %f, %f}\n", y[0],y[1],y[2]);
 }

 This will replace y by: {-2.0, 2.0, 0.0} and return meam=3.0.

 In the particular call:

 removeVectorMean(myModel.AccessRPTemplate()+j*nv_myModel.AccessRPTemplate(), nv_myModel.AccessRPTemplate());

 The array y is represented by: myModel.AccessRPTemplate()+j*nv_myModel.AccessRPTemplate(), and n =
 nv_myModel.AccessRPTemplate()
 */
template <typename ValuesType>
ValuesType
removeVectorMean(std::vector<ValuesType> & x)
{
  const ValuesType n = x.size();
  const ValuesType mean = std::accumulate(x.begin(), x.end(), static_cast<ValuesType>(0.0)) / n;
  for (auto it = x.begin(); it != x.end(); ++it)
  {
    *it -= mean;
  }
  return mean;
}

/**
 * This function takes an image and returns an array of type
 float comprised of an ordered subset of the intensity values
 stored in 'image'.  The subset of  voxels that is returned is
 centered around point 'CenterPoint'.  The shape of the subset is
 specified by a set of point offesets from the 'CenterPoint'
 given specified by the vector 'model'.

 If any computed voxel location happens to fall outside the
 memory block 'image', a value of 0.0 is returned in for for that
 location in 'result_array'.
 *
 * @author hjohnson (8/29/2008)
 *
 * @param image
 * @param CenterPoint
 * @param model
 * @param result_array
 */
extern void
extractArray(const LinearInterpolatorType::Pointer &                        imInterp,
             const SImageType::PointType &                                  CenterPoint,
             const landmarksConstellationModelIO::IndexLocationVectorType & model,
             std::vector<float> &                                           result_array);

#include <itkMinimumMaximumImageFilter.h>
#include <itkScalarImageToHistogramGenerator.h>
#include <itkOtsuMultipleThresholdsCalculator.h>
template <typename SImageType>
void
ImageMinMax(typename SImageType::Pointer     image,
            typename SImageType::PixelType * imageMin,
            typename SImageType::PixelType * imageMax)
{
  typename itk::MinimumMaximumImageFilter<SImageType>::Pointer minmaxFilter =
    itk::MinimumMaximumImageFilter<SImageType>::New();
  minmaxFilter->SetInput(image);
  minmaxFilter->Update();
  *imageMax = minmaxFilter->GetMaximum();
  *imageMin = minmaxFilter->GetMinimum();
}

/**
 *
 *
 * @author hjohnson (8/12/2008)
 *
 * @param SImageType The image type templated over
 * @param image  the itk image to be used to compute the
 *               histograms
 * @param low   The intensity value where "percent" voxels are
 *              above this threshold
 * @param high  The intensity value where "1.0-percent" voxels
 *              are below this threshold
 * @param percent A value between 0.0 and 100.0 representing the
 *                Quantile percentages to be eliminated. NOTE:
 *                This value will be divided by 100.0, so
 *                100.0=100%, and 1.0=.01%.
 * @return Background threshold.
 */
template <typename SImageType>
typename SImageType::PixelType
setLowHigh(typename SImageType::Pointer &   image,
           typename SImageType::PixelType & low,
           typename SImageType::PixelType & high,
           const float                      percent)
{
  typename SImageType::PixelType imageMin;
  typename SImageType::PixelType imageMax;
  ImageMinMax<SImageType>(image, &imageMin, &imageMax);

  using HistogramGeneratorType = itk::Statistics::ScalarImageToHistogramGenerator<SImageType>;
  typename HistogramGeneratorType::Pointer histogramGenerator = HistogramGeneratorType::New();
  histogramGenerator->SetInput(image);
  histogramGenerator->SetNumberOfBins(imageMax - imageMin + 1);
  histogramGenerator->SetMarginalScale(1.0);
  histogramGenerator->SetHistogramMin(imageMin);
  histogramGenerator->SetHistogramMax(imageMax);
  histogramGenerator->Compute();
  using HistogramType = typename HistogramGeneratorType::HistogramType;
  const HistogramType * histogram = histogramGenerator->GetOutput();

  using OtsuCalcType = itk::OtsuMultipleThresholdsCalculator<HistogramType>;
  typename OtsuCalcType::Pointer OtsuCalc = OtsuCalcType::New();
  OtsuCalc->SetInputHistogram(histogram);
  OtsuCalc->SetNumberOfThresholds(1);
  OtsuCalc->Compute();
  typename OtsuCalcType::OutputType otsuThresholds = OtsuCalc->GetOutput();

  low = static_cast<typename SImageType::PixelType>(histogram->Quantile(0, 0.0F + percent));
  high = static_cast<typename SImageType::PixelType>(histogram->Quantile(0, 1.0F - percent));
  return static_cast<typename SImageType::PixelType>(otsuThresholds[0]);
}




extern SimilarityTransformType::Pointer
DoIt_Similarity(const PointList & fixedPoints, const PointList & movingPoints);

extern VersorRigidTransformType::Pointer
ComputeRigidTransformFromLandmarkLists(const PointList & fixedPoints, const PointList & movingPoints);

#endif //__landmarksConstellationCommon_h
