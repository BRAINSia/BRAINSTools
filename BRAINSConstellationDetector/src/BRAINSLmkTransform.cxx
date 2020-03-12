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
 * Author: Wei Lu
 * at Psychiatry Imaging Lab, University of Iowa Health Care, 2010
 */

#include "itkImage.h"
#include "itkPoint.h"
#include "itkPointSet.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkResampleImageFilter.h"
#include "itkAffineTransform.h"
#include "itkTransformFileWriter.h"
#include "itkLinearInterpolateImageFunction.h"
#include "BRAINSThreadControl.h"

// Use modified itkKernelTransform to get affine transform
#include "itkKernelTransform.h"
#include "itkThinPlateSplineKernelTransform.h"
#include <BRAINSCommonLib.h>

template <typename TScalarType, unsigned int NDimension>
class BCDThinPlateSplineKernelTransform : public itk::ThinPlateSplineKernelTransform<TScalarType, NDimension>
{
public:
  /** Standard class type alias. */
  using Self = BCDThinPlateSplineKernelTransform;
  using Superclass = itk::ThinPlateSplineKernelTransform<TScalarType, NDimension>;
  using Pointer = itk::SmartPointer<Self>;
  using ConstPointer = itk::SmartPointer<const Self>;

  /** New macro for creation of through a Smart Pointer */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(BCDThinPlateSplineKernelTransform, ThinPlateSplineKernelTransform);
  typename Superclass::Superclass::AMatrixType
  GetAMatrix()
  {
    return this->Superclass::Superclass::m_AMatrix;
  }

  typename Superclass::Superclass::BMatrixType
  GetBVector()
  {
    return this->Superclass::Superclass::m_BVector;
  }
};

#include "BRAINSLmkTransformCLP.h"

#include <fstream>
#include <vector>
#include <iostream>

/*
 * Description:
 *
 * This utility program estimates the affine transform to align the fixed landmarks
 * to the moving landmarks, and then generate the resampled moving image to the same
 * physical space as that of the reference image
 */

constexpr unsigned int ImageDimension = 3;
using PixelType = short;
using ImageType = itk::Image<PixelType, ImageDimension>;
using LandmarksVectorType = std::vector<ImageType::PointType>;

LandmarksVectorType
LoadLandmarks(std::string filename);

int
main(int argc, char * argv[])
{
  PARSE_ARGS;
  BRAINSRegisterAlternateIO();

  const BRAINSUtils::StackPushITKDefaultNumberOfThreads TempDefaultNumberOfThreadsHolder(numberOfThreads);
  if ((inputMovingLandmarks.compare("") == 0) && (inputFixedLandmarks.compare("") == 0) &&
      (inputMovingVolume.compare("") == 0) && (inputReferenceVolume.compare("") == 0))
  {
    itkGenericExceptionMacro(<< "Please set inputMovingLandmarks, inputFixedLandmarks, "
                             << "inputMovingVolume, and inputReferenceVolume.");
  }

  // type alias
  using CoordinateRepType = double;
  using ImageReaderType = itk::ImageFileReader<ImageType>;
  using ImageWriterType = itk::ImageFileWriter<ImageType>;
  using TPSTransformType = BCDThinPlateSplineKernelTransform<CoordinateRepType, ImageDimension>;
  using AffineTransformType = itk::AffineTransform<CoordinateRepType, ImageDimension>;
  using PointSetType = TPSTransformType::PointSetType;
  using TransformWriterType = itk::TransformFileWriter;
  using PointIdType = PointSetType::PointIdentifier;
  using ResamplerType = itk::ResampleImageFilter<ImageType, ImageType>;
  using InterpolatorType = itk::LinearInterpolateImageFunction<ImageType, double>;

  // Read in landmarks
  PointSetType::Pointer sourceLandmarks = PointSetType::New();
  PointSetType::Pointer targetLandmarks = PointSetType::New();
  {
    PointSetType::PointsContainer::Pointer sourceLandmarkContainer = sourceLandmarks->GetPoints();
    PointSetType::PointsContainer::Pointer targetLandmarkContainer = targetLandmarks->GetPoints();
    PointIdType                            id = itk::NumericTraits<PointIdType>::ZeroValue();
    LandmarksVectorType                    targetLandmarksVec = LoadLandmarks(inputMovingLandmarks);
    LandmarksVectorType                    sourceLandmarksVec = LoadLandmarks(inputFixedLandmarks);

    // Sanity check
    if (targetLandmarksVec.size() != sourceLandmarksVec.size())
    {
      std::cerr << "Different number of fixed and moving landmarks!" << std::endl;
      return EXIT_FAILURE;
    }

    unsigned int idx = 0;
    for (idx = 0; idx < sourceLandmarksVec.size(); ++idx)
    {
      sourceLandmarkContainer->InsertElement(id, sourceLandmarksVec[idx]);
      targetLandmarkContainer->InsertElement(id++, targetLandmarksVec[idx]);
    }
  }

  // Estimate affine transform
  AffineTransformType::Pointer affine = AffineTransformType::New();
  {
    TPSTransformType::Pointer tps = TPSTransformType::New();
    tps->SetSourceLandmarks(sourceLandmarks);
    tps->SetTargetLandmarks(targetLandmarks);
    tps->ComputeWMatrix();
    itk::Matrix<double, ImageDimension, ImageDimension> aMatrix(tps->GetAMatrix());
    itk::Vector<double, ImageDimension>                 bVector;
    bVector.SetVnlVector(vnl_vector<double>(tps->GetBVector()));
    itk::Matrix<double, ImageDimension, ImageDimension> identity;
    identity.SetIdentity();
    affine->SetMatrix(aMatrix + identity);
    affine->SetOffset(bVector);
  }

  // Write output aligning transform
  if (outputAffineTransform.compare("") != 0)
  {
    TransformWriterType::Pointer writer = TransformWriterType::New();
    writer->SetInput(affine);
    writer->SetFileName(outputAffineTransform);
    writer->SetUseCompression(true);
    try
    {
      writer->Update();
    }
    catch (itk::ExceptionObject & excep)
    {
      std::cerr << "Cannot write the outputTransform file!" << std::endl;
      std::cerr << excep << std::endl;
    }
  }

  // Read in images
  ImageType::Pointer movingImage;
  {
    ImageReaderType::Pointer reader = ImageReaderType::New();
    reader->SetFileName(inputMovingVolume);
    try
    {
      reader->Update();
    }
    catch (itk::ExceptionObject & excp)
    {
      std::cerr << "Exception thrown " << std::endl;
      std::cerr << excp << std::endl;
      return EXIT_FAILURE;
    }
    movingImage = reader->GetOutput();
  }

  ImageType::Pointer referenceImage;
  {
    ImageReaderType::Pointer reader = ImageReaderType::New();
    reader->SetFileName(inputReferenceVolume);
    try
    {
      reader->Update();
    }
    catch (itk::ExceptionObject & excp)
    {
      std::cerr << "Exception thrown " << std::endl;
      std::cerr << excp << std::endl;
      return EXIT_FAILURE;
    }
    referenceImage = reader->GetOutput();
  }

  // Resample moving image
  ResamplerType::Pointer resampler = ResamplerType::New();
  {
    InterpolatorType::Pointer interpolator = InterpolatorType::New();
    resampler->SetUseReferenceImage(true);
    resampler->SetInput(movingImage);
    resampler->SetReferenceImage(referenceImage);
    resampler->SetInterpolator(interpolator);
    resampler->SetTransform(affine);
  }

  // Write aligned image
  if (outputResampledVolume.compare("") != 0)
  {
    ImageWriterType::Pointer writer = ImageWriterType::New();
    writer->SetInput(resampler->GetOutput());
    writer->SetFileName(outputResampledVolume);
    writer->SetUseCompression(true);
    try
    {
      writer->Update();
    }
    catch (itk::ExceptionObject & excp)
    {
      std::cerr << "Exception thrown " << std::endl;
      std::cerr << excp << std::endl;
      return EXIT_FAILURE;
    }
  }

  return EXIT_SUCCESS;
}

LandmarksVectorType
LoadLandmarks(std::string filename)
{
  LandmarksVectorType landmarks;
  std::string         line;
  std::ifstream       myfile(filename.c_str());

  if (!myfile.is_open())
  {
    itkGenericExceptionMacro(<< "Fatal error: Failed to load landmarks file. Program abort!");
  }
  while (getline(myfile, line))
  {
    if (line.compare(0, 1, "#") != 0)
    {
      unsigned int i = 0;
      int          pos1 = line.find(',', 0);
      int          pos2 = 0;
      std::string  name = line.substr(0, pos1);
      if (name.compare("CM") == 0) // exclude CM
      {
        continue;
      }
      ImageType::PointType labelPos;
      for (i = 0; i < 3; ++i)
      {
        pos2 = line.find(',', pos1 + 1);
        labelPos[i] = std::stod(line.substr(pos1 + 1, pos2 - pos1 - 1).c_str());
        if (i < 2)
        {
          labelPos[i] *= -1; // RAS -> LPS
        }
        pos1 = pos2;
      }
      landmarks.push_back(labelPos);
    }
  }

  myfile.close();
  return landmarks;
}
