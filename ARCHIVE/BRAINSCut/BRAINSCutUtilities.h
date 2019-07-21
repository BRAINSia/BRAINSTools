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
#ifndef BRAINSCutUtilities_h
#define BRAINSCutUtilities_h

#include "BRAINSCutConfiguration.h"

#include "TrainingVectorConfigurationType.h"
#include "TrainingPrameters.h"
#include "ApplyModel.h"
#include "itkIO.h"

#include "GenericTransformImage.h"

#include <itkSmoothingRecursiveGaussianImageFilter.h>

/** include opencv library */
#include "opencv2/ml.hpp"

// #include "opencv2/flann/flann.hpp" **opencv 2
#include "opencv2/flann.hpp"

// using OpenCVMLPType = CvANN_MLP_Revision;
// using OpenCVMLPType = CvANN_MLP;
using OpenCVMLPType = cv::ml::ANN_MLP;

/** Training data set definition */

using scalarType = float;

struct pairedTrainingSetType
{
  cv::Mat      pairedInput;
  cv::Mat      pairedOutput;
  cv::Mat      pairedOutputRF;
  unsigned int size;
};

/* normalization type */
enum FeatureNormalizationMethodEnum
{
  Linear,
  Sigmoid,
  DoubleSigmoid,
  zScore
};
/*
 * constant
 */
static const float            HundredPercentValue = 1.0F;
static const float            ZeroPercentValue = -1.0F;
static constexpr unsigned int LineGuardSize = 1;
static const scalarType       LineGuard = 1234567.0;
static const float            FLOAT_TOLERANCE = 0.01;

/*
 * Image Definitions
 */
constexpr unsigned char DIMENSION = 3;

using ReadInPixelType = double;
using ReadInImageType = itk::Image< ReadInPixelType, DIMENSION >;
using ReadInImagePointer = ReadInImageType::Pointer;

using WorkingPixelType = float;
using WorkingImageType = itk::Image< WorkingPixelType, DIMENSION >;
using WorkingImagePointer = WorkingImageType::Pointer;

using WorkingImageVectorType = std::vector< WorkingImagePointer >;

/* Deformations */
using DeformationScalarType = double;
using DeformationPixelType = itk::Vector< DeformationScalarType, DIMENSION >;
using DisplacementFieldType = itk::Image< DeformationPixelType, DIMENSION >;

using WorkingIndexType = WorkingImageType::IndexType;

using InputVectorType = std::vector< WorkingPixelType >;
using OutputVectorType = std::vector< WorkingPixelType >;

// HACK INFO:  Regina int below should be unsigned int to avoid negative index numbers
using InputVectorMapType = std::map< int, InputVectorType >; // < index ,feature vector > pair
using OutputVectorMapType = std::map< int, OutputVectorType >;
using PredictValueMapType = std::map< int, scalarType >;

std::string
GetAtlasToSubjectRegistrationFilename( DataSet & subject );

std::string
GetSubjectToAtlasRegistrationFilename( DataSet & subject );

WorkingImagePointer
SmoothImage( const WorkingImagePointer & image, const float GaussianValue );

WorkingImagePointer
ReadImageByFilename( const std::string & filename );

DisplacementFieldType::Pointer
GetDeformationField( std::string filename );

itk::Transform< double, 3, 3 >::Pointer
GetGenericTransform( std::string filename );

#endif
