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
#include "ml.h"
#include "cxcore.h"

// #include "opencv2/flann/flann.hpp" **opencv 2
#include "opencv2/flann.hpp"

// typedef CvANN_MLP_Revision OpenCVMLPType;
//typedef CvANN_MLP OpenCVMLPType;
typedef cv::ml::ANN_MLP OpenCVMLPType;

/** Training data set definition */

typedef float scalarType;

struct pairedTrainingSetType
  {
  cv::Mat pairedInput;
  cv::Mat pairedOutput;
  cv::Mat pairedOutputRF;
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
static const float        HundredPercentValue = 1.0F;
static const float        ZeroPercentValue = -1.0F;
static constexpr unsigned int LineGuardSize = 1;
static const scalarType   LineGuard = 1234567.0;
static const float        FLOAT_TOLERANCE = 0.01;

/*
* Image Definitions
*/
constexpr unsigned char DIMENSION = 3;

typedef double                                 ReadInPixelType;
typedef itk::Image<ReadInPixelType, DIMENSION> ReadInImageType;
typedef ReadInImageType::Pointer               ReadInImagePointer;

typedef float                                   WorkingPixelType;
typedef itk::Image<WorkingPixelType, DIMENSION> WorkingImageType;
typedef WorkingImageType::Pointer               WorkingImagePointer;

typedef std::vector<WorkingImagePointer> WorkingImageVectorType;

/* Deformations */
typedef double                                        DeformationScalarType;
typedef itk::Vector<DeformationScalarType, DIMENSION> DeformationPixelType;
typedef itk::Image<DeformationPixelType, DIMENSION>   DisplacementFieldType;

typedef WorkingImageType::IndexType WorkingIndexType;

typedef std::vector<WorkingPixelType> InputVectorType;
typedef std::vector<WorkingPixelType> OutputVectorType;

// HACK TODO:  Regina int below should be unsigned int to avoid negative index numbers
typedef std::map<int, InputVectorType>  InputVectorMapType; // < index ,feature vector > pair
typedef std::map<int, OutputVectorType> OutputVectorMapType;
typedef std::map<int, scalarType>       PredictValueMapType;

std::string GetAtlasToSubjectRegistrationFilename( DataSet& subject);

std::string GetSubjectToAtlasRegistrationFilename( DataSet& subject);

WorkingImagePointer SmoothImage( const WorkingImagePointer& image, const float GaussianValue);

WorkingImagePointer ReadImageByFilename( const std::string  & filename );

DisplacementFieldType::Pointer GetDeformationField( std::string filename);

itk::Transform<double, 3, 3>::Pointer GetGenericTransform( std::string filename);

#endif
