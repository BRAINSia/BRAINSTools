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
 * This is a temporary file to be used during refactoring to help consolicate a family of functions into a common library.
 */
#ifndef __landmarksConstellationCommon_h
#define __landmarksConstellationCommon_h

// Use linear interpolation to keep the processing quick.
//RM #define __QUICK_RUNS_APPROXIMATION__

#include <cstdio>      // TODO: This include file should be removed, prefer constructs
                       // from the std library
#include <cstdlib>     // TODO: This include file should be removed, prefer constructs
                       // from the std library
#include <cmath>       // TODO: This include file should be removed, use vcl_math
                       // instead
#include <sys/types.h> // TODO: This include file should be removed, unix only
                       // non-portable
#include <sys/stat.h>  // TODO: This include file should be removed, unix only
                       // non-portable
#include <unistd.h>    // TODO: This include file should be removed, unix only
                       // non-portable
#include <ctime>       // TODO: This include file should be removed, unix only
                       // non-portable
#include <cctype>      // TODO: This include file should be removed, use vcl_math
                       // instead
                       // #include <volume.h> //This include file should be
                       // removed

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

#include "itkCylinderSpatialObject.h"
#include "itkSpatialObjectToImageFilter.h"
#include "itkRGBPixel.h"
#include "itkIO.h"
#include "itkExtractImageFilter.h"
#include "itkStatisticsImageFilter.h"
#include "itkFlipImageFilter.h"
#include <vcl_compiler.h>
#include <iostream>
#include "algorithm"
#include "itkImageDuplicator.h"
#include "itkEuler3DTransform.h"
#include <vnl/vnl_cross.h>

#include <itkScalarImageToHistogramGenerator.h>

#include "GenericTransformImage.h"

#include "Slicer3LandmarkIO.h"

#if 0 //RM
extern const unsigned int MAX_ROTATIONS_TESTED;
extern const unsigned int MAXITER;
extern const unsigned int DEL;
extern const unsigned int YES;
extern const unsigned int SMAX;
#endif
extern const unsigned int NO;
namespace LMC
{
extern bool debug;
extern bool globalverboseFlag;
}

//
//
// ////////////////////////////////////////////////////////////////////////////////////////////////

namespace // avoid 'shadows declaration' warnings.
{
//RM typedef float                 vertexType[4][3];
typedef itk::Image<short, 3>  SImageType;
typedef itk::Image<double, 3> DImageType3D;
typedef itk::Image<float, 3>  FImageType3D;
//RM typedef itk::Image<short,2> SImageType2D;
//RM typedef itk::Image<double, 2> DImageType2D;
//RM typedef itk::Image<float, 2>  FImageType2D;

typedef itk::Image<unsigned char, 3> ByteImageType;

typedef SImageType::PointType                   SImagePointType;

typedef itk::RGBPixel<unsigned char> RGBPixelType;
typedef itk::Image<RGBPixelType, 3>  RGBImageType;
typedef itk::Image<RGBPixelType, 2>  RGB2DImageType;

typedef itk::Euler3DTransform<double>       RigidTransformType;
typedef itk::VersorRigid3DTransform<double> VersorTransformType;

typedef itk::MultiResolutionPyramidImageFilter<SImageType, SImageType> PyramidFilterType;
typedef itk::LinearInterpolateImageFunction<SImageType, double>        LinearInterpolatorType;
}

#include "landmarksConstellationModelIO.h"

//RM extern VersorTransformType::Pointer ConvertToVersorRigid3D(RigidTransformType::Pointer RT);

extern std::string globalResultsDir;
extern int         globalImagedebugLevel;
extern PyramidFilterType::Pointer MakeThreeLevelPyramid(SImageType::Pointer refImage);

extern PyramidFilterType::Pointer MakeOneLevelPyramid(SImageType::Pointer refImage);

extern SImageType::PointType GetImageCenterPhysicalPoint(SImageType::Pointer & image);

extern unsigned char ShortToUChar(short in, short min, short max);

extern SImageType::Pointer CreateTestCenteredRotatedImage2(const RigidTransformType::Pointer ACPC_MSP_AlignedTransform,
                                                           /* const
                                                             SImageType::PointType
                                                             finalPoint, */
                                                           const SImageType::PointType PreMSP_Point,
                                                           /*const*/ SImageType::Pointer & image,
                                                           const RigidTransformType::Pointer & Point_Rotate);

#if 0 //RM
extern itk::Matrix<double, 3, 3> GetMatrixInverse(const itk::Matrix<double, 3, 3> & input);

// extern itk::Matrix<double,3,3> CreateRotationMatrixFromAngles(const double
// alpha, const double beta, const double gamma);
extern itk::Versor<double> CreateRotationVersorFromAngles(const double alpha, const double beta, const double gamma);

extern void ComputeEulerAnglesFromRotationMatrix(const itk::Matrix<double, 3,
                                                                   3> &  m, double & initialAttitudeAngle,
                                                 double & initialBankAngle, double & initialHeadingAngle);

extern void defineTemplateIndexLocations(const int r, const int h,
                                         landmarksConstellationModelIO::IndexLocationVectorType & indexLocations,
                                         std::vector<float> & results_array);

extern int computeTemplateSize(const int r, const int h);

#endif
extern void decomposeRPAC(const SImageType::PointType & RP, const SImageType::PointType & PC,
                          const SImageType::PointType & AC, double *const RPPC_to_RPAC_angle,
                          double *const RPAC_over_RPPC);

extern void MakeLabelImage(SImageType::Pointer in, const SImageType::PointType & RP, const SImageType::PointType & AC,
                           const SImageType::PointType & PC, const SImageType::PointType & VN4,
                           const std::string & fname);

#if 0 //RM
extern SImageType::PointType::VectorType initialAC(const SImageType::PointType & RP, const SImageType::PointType & PC,
                                                   const double RPPC_to_RPAC_angleMean,
                                                   const double RPAC_over_RPPCMean);
#endif

//RM typedef itk::Statistics::MersenneTwisterRandomVariateGenerator RandomGeneratorType;
//RM typedef RandomGeneratorType::Pointer                           RandomGeneratorPointer;
//RM extern RandomGeneratorPointer _RandomGenerator;
//RM extern double GetRandomZeroOneDouble(void);

//RM extern void InitializeRandomZeroOneDouble(RandomGeneratorType::IntegerType rseed);

extern void ComputeMSP(SImageType::Pointer image, RigidTransformType::Pointer & Tmsp,
                       SImageType::Pointer & transformedImage, const SImageType::PointType & centerOfHeadMass,
                       const int qualityLevel, double & cc);

extern void ComputeMSP_Easy(SImageType::Pointer image, RigidTransformType::Pointer & Tmsp,
                            const SImageType::PointType & centerOfHeadMass, const int qualityLevel);

extern SImageType::Pointer CreatedebugPlaneImage(SImageType::Pointer referenceImage,
                                                 const RigidTransformType::Pointer MSPTransform,
                                                 const std::string & debugfilename);

extern void CreatedebugPlaneImage(SImageType::Pointer referenceImage, const std::string & debugfilename);

extern RigidTransformType::Pointer computeTmspFromPoints(SImageType::PointType RP, SImageType::PointType AC,
                                                         SImageType::PointType PC, SImageType::PointType DesiredCenter);

extern SImageType::PointType GetCenterOfHeadMass(SImageType::Pointer volume);

extern SImageType::Pointer MakeIsoTropicReferenceImage();

/*****************************************************************************/
template <class ValuesType>
ValuesType vectorNorm(const std::vector<ValuesType> & x)
{
  ValuesType norm = 0.0;

  for( typename std::vector<ValuesType>::const_iterator it = x.begin();
       it != x.end(); ++it )
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
template <class ValuesType>
void normalizeVector(std::vector<ValuesType> & x)
{
  const ValuesType norm = vectorNorm(x);

  if( norm < std::numeric_limits<ValuesType>::epsilon() )
    {
    std::cout << "WARNING:  ZERO NORM VECTOR." << __FILE__ << __LINE__ << std::endl;
    return;
    }
  for( typename std::vector<ValuesType>::iterator it = x.begin();
       it != x.end(); ++it )
    {
    *it /= norm;
    }
  return;
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
template <class ValuesType>
ValuesType removeVectorMean(std::vector<ValuesType> & x)
{
  const ValuesType n = x.size();
  const ValuesType mean = std::accumulate(x.begin(),x.end(),static_cast<ValuesType>(0.0) ) / n;
  for( typename std::vector<ValuesType>::iterator it = x.begin();
       it != x.end(); ++it )
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
extern
void extractArray(
  LinearInterpolatorType::Pointer imInterp,
  const SImageType::PointType & CenterPoint,
  const landmarksConstellationModelIO::IndexLocationVectorType & model,
  std::vector<float> & result_array);

#if 0 //RM
inline
static std::string
PrefixName(const char *prefix, const std::string & name)
{
  std::string rval;

  rval += itksys::SystemTools::GetFilenamePath(name);
  if( rval.size() > 0 )
    {
    rval += "/";
    }
  //  std::string rval(pathpart);
  rval += prefix;
  rval += itksys::SystemTools::GetFilenameName(name);
  return rval;
}
#endif

#include <itkMinimumMaximumImageFilter.h>
#include <itkScalarImageToHistogramGenerator.h>
#include <itkOtsuMultipleThresholdsCalculator.h>
template <class SImageType>
void ImageMinMax(typename SImageType::Pointer image,
                 typename SImageType::PixelType *imageMin, typename SImageType::PixelType *imageMax)
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
template <class SImageType>
typename SImageType::PixelType
setLowHigh(typename SImageType::Pointer & image,
           typename SImageType::PixelType & low,
           typename SImageType::PixelType & high,
           const float percent)
{
  typename SImageType::PixelType imageMin;
  typename SImageType::PixelType imageMax;
  ImageMinMax<SImageType>(image, &imageMin, &imageMax);

  typedef itk::Statistics::ScalarImageToHistogramGenerator<SImageType> HistogramGeneratorType;
  typename HistogramGeneratorType::Pointer histogramGenerator = HistogramGeneratorType::New();
  histogramGenerator->SetInput(image);
  histogramGenerator->SetNumberOfBins(imageMax - imageMin + 1);
  histogramGenerator->SetMarginalScale(1.0);
  histogramGenerator->SetHistogramMin(imageMin);
  histogramGenerator->SetHistogramMax(imageMax);
  histogramGenerator->Compute();
  typedef typename HistogramGeneratorType::HistogramType HistogramType;
  const HistogramType *histogram = histogramGenerator->GetOutput();

  typedef itk::OtsuMultipleThresholdsCalculator<HistogramType> OtsuCalcType;
  typename OtsuCalcType::Pointer OtsuCalc = OtsuCalcType::New();
  OtsuCalc->SetInputHistogram(histogram);
  OtsuCalc->SetNumberOfThresholds(1);
  OtsuCalc->Update();
  typename OtsuCalcType::OutputType otsuThresholds = OtsuCalc->GetOutput();

  low = static_cast<typename SImageType::PixelType>( histogram->Quantile(0, 0.0F + percent) );
  high = static_cast<typename SImageType::PixelType>( histogram->Quantile(0, 1.0F - percent) );
  return static_cast<typename SImageType::PixelType>( otsuThresholds[0] );
}

#if 0 //TODO: Remove old code
// ------------------------------
// The following should be cleaned up and moved elsewhere
template <class DType>
double
removeVectorMean(DType *x, DType *y, int n)
{
  double mean = 0.0;

  if( n <= 0 )
    {
    return 0.0;
    }
  for( int i = 0; i < n; ++i )
    {
    mean += static_cast<double>( x[i] );
    }
  mean /= n;
  for( int i = 0; i < n; ++i )
    {
    y[i] = static_cast<DType>( x[i] - mean );
    }
  return mean;
}

extern void removeVectorMean(double *y, int n, int p);

template <class DType>
double removeVectorMean(DType *y, int n)
{
  double mean = 0.0;

  if( n <= 0 )
    {
    return 0.0;
    }
  for( int i = 0; i < n; ++i )
    {
    mean += y[i];
    }
  mean /= n;
  for( int i = 0; i < n; ++i )
    {
    y[i] = y[i] - mean;
    }
  return mean;
}

// ///////////////////////////////////////////////////////
// Implements Eq. (2.2.1) of J. Cohen & P. Cohen (2nd ed.)
template <class DType>
double standardDeviation(DType *x, int n)
{
  double sx, sxx;
  double sd;

  if( n <= 0 )
    {
    return 0.0;
    }
  sx = sxx = 0.0;
  for( int i = 0; i < n; ++i )
    {
    sx += x[i];
    sxx += x[i] * x[i];
    }
  sd = ( sxx - sx * sx / n ) / n;
  if( sd > 0.0 )
    {
    sd = std::sqrt(sd);
    }
  else
    {
    sd = 0.0;
    }
  return sd;
}

// Computes the Pearson correlation coefficient between x and y.
// Implements Eq. (2.3.2) of J. Cohen & P. Cohen (2nd ed.)
template <class DTypeX, class DTypeY>
double pearsonCorrelation(DTypeX *x, DTypeY *y, int n)
{
  double sx, sy, sxx, syy, sxy;
  double Sxx, Sxy, Syy;
  double dum = 0.0;

  sx = sy = sxx = syy = sxy = 0.0;
  for( int i = 0; i < n; ++i )
    {
    sx += x[i];
    sxx += x[i] * x[i];
    sy += y[i];
    syy += y[i] * y[i];
    sxy += x[i] * y[i];
    }
  Sxx = n * sxx - sx * sx;
  Syy = n * syy - sy * sy;
  Sxy = n * sxy - sx * sy;
  if( Sxx * Syy > 0.0 )
    {
    dum = std::sqrt(Sxx * Syy);
    }
  if( dum != 0.0 )
    {
    return (double)( Sxy / dum );
    }
  else
    {
    return 0.0;
    }
}

// Implements Eq. (3.3.11) of J. Cohen & P. Cohen (2nd ed.)
template <class DTypeY, class DTypeX>
void partialCorrelation(DTypeY *Y, DTypeX *X1, DTypeX *X2, int n, double *pr1, double *pr2)
{
  double rY1, rY2, r12;
  double dum1, dum2;

  rY1 = pearsonCorrelation(X1, Y, n);
  rY2 = pearsonCorrelation(X2, Y, n);
  r12 = pearsonCorrelation(X1, X2, n);
  dum1 = ( 1.0 - rY2 * rY2 ) * ( 1 - r12 * r12 );
  dum2 = ( 1.0 - rY1 * rY1 ) * ( 1 - r12 * r12 );
  if( dum1 > 0.0 )
    {
    *pr1 = ( rY1 - rY2 * r12 ) / std::sqrt(dum1);
    }
  else
    {
    *pr1 = 0.0;
    }
  if( dum2 > 0.0 )
    {
    *pr2 = ( rY2 - rY1 * r12 ) / std::sqrt(dum2);
    }
  else
    {
    *pr2 = 0.0;
    }
}

// //////////////////////////////////////////////////////////////
// Computes the sample mean of a set of n observations {x_1,x_2,...,x_n} from a
// given distribution.
// //////////////////////////////////////////////////////////////
template <typename DType>
double
sample_mean(DType *x, int n)
{
  double mean = 0.0;

  if( n <= 0 )
    {
    return 0.0;
    }
  for( int i = 0; i < n; ++i )
    {
    mean += x[i];
    }
  mean /= n;
  return mean;
}

// //////////////////////////////////////////////////////////////
// Computes the unbiased sample variance of a set of n observations
// {x_1,x_2,...,x_n} from a given distribution.
// //////////////////////////////////////////////////////////////
template <typename DType>
double sample_variance(DType *x, int n, double *mean)
{
  double sum_of_sq = 0.0;
  double sum = 0.0;
  double var;

  if( n < 2 )
    {
    *mean = 0.0;
    return 0.0;
    }
  for( int i = 0; i < n; ++i )
    {
    sum_of_sq +=
      static_cast<double>( x[i] )
      * static_cast<double>( x[i] );
    sum += static_cast<double>( x[i] );
    }
  *mean = sum / n;
  var = ( sum_of_sq - sum * sum / n ) / ( n - 1.0 );
  return var;
}

template <class DType>
double
independent_samples_t(DType *x1, int n1, DType *x2, int n2, int *df, double *meandiff)
{
  double var1, var2;
  double mean1, mean2;
  double var, sd;
  double t;

  *df = n1 + n2 - 2;
  if( n1 < 2 || n2 < 2 )
    {
    return 0.0;
    }
  var1 = sample_variance<DType>(x1, n1, &mean1);
  var2 = sample_variance<DType>(x2, n2, &mean2);
  var = ( ( ( n1 - 1 ) * var1 + ( n2 - 1 ) * var2 ) / ( n1 + n2 - 2.0 ) ) * ( 1.0 / n1 + 1.0 / n2 );
  sd = std::sqrt(var);
  if( sd == 0.0 )
    {
    return 0.0;
    }
  *meandiff = mean1 - mean2;
  t = ( mean1 - mean2 ) / sd;
  return t;
}

template <class DType>
double
paired_samples_t(DType *x1, DType *x2, int n, int *df, double *meandiff)
{
  double var;
  double mean;
  double sd;
  double t;

  *df = n - 1;
  if( n < 2 )
    {
    return 0.0;
    }
  for( int i = 0; i < n; ++i )
    {
    x1[i] -= x2[i];
    }
  var = sample_variance<double>(x1, n, &mean) / n;
  sd = std::sqrt(var);
  if( sd == 0.0 )
    {
    return 0.0;
    }
  *meandiff = mean;
  t = mean / sd;
  return t;
}
#endif
#if 0
template<class TScalarType>
extern void WriteTransformToDisk( itk::Transform<TScalarType, 3, 3> * myTransform , const std::string & filename  );
#endif

#endif
