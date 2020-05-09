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
#include "itkLandmarkBasedTransformInitializer.h"
#include "itkVersorRigid3DTransform.h"
#include "itkSimilarity3DTransform.h"
#include "itkAffineTransform.h"
#include "itkBSplineTransform.h"
#include "itkImage.h"
#include "itkImageFileReader.h"

#include <map>

#include "GenericTransformImage.h"
#include "Slicer3LandmarkIO.h"

#include "BRAINSLandmarkInitializerCLP.h"

using ParameterValueType = double;
using PixelType = double;
constexpr unsigned int Dimension = 3;
using ImageType = itk::Image<PixelType, Dimension>;

using LandmarkWeightType = std::vector<double>;
using LandmarkWeightConstIterator = LandmarkWeightType::const_iterator;
using LandmarkPointType = itk::Point<double, Dimension>;
using LandmarkPointContainer = std::vector<LandmarkPointType>;


static void
CheckLandmarks(const LandmarksMapType & ldmk, const LandmarksWeightMapType & weightMap)
{
  if (ldmk.size() < 3)
  {
    std::cerr << "At least 3 fiducial points must be specified. " << std::endl;
    exit(EXIT_FAILURE);
  }

  //  if( ldmk.find( "AC" ) == ldmk.end() ||
  //      ldmk.find( "PC" ) == ldmk.end() ||
  //      ldmk.find( "LE" ) == ldmk.end() ||
  //      ldmk.find( "RE" ) == ldmk.end() )
  //    {
  //    std::cerr << " Base four landmarks ( AC, PC, left eye(LE), and right eye(RE) ) "
  //              << " has to be provided"
  //              << std::endl;
  //    exit(EXIT_FAILURE);
  //    }

  for (const auto & i : weightMap)
  {
    if (ldmk.find(i.first) == ldmk.end())
    {
      std::cerr << "WARNING: Landmark not found: " << i.first << std::endl;
    }
#if defined(VERBOSE_OUTPUT)
    else
    {
      std::cerr << "NOTE: Landmark found: " << i->first << std::endl;
    }
#endif
  }
}

double ComputeIsotropicScaleFactor(
  const LandmarkPointContainer & fixedLmks,
  const LandmarkPointContainer & movingLmks,
  const LandmarkWeightType & landmarkWgts
  )
{
  double isotropic_scale_factor=1.0;
  const size_t numLmks =fixedLmks.size();

  const bool has_lmk_weights = ( landmarkWgts.size() == numLmks);

  double fix_dist_sum = 0.0;
  double mov_dist_sum = 0.0;
  const auto & fix_ref = fixedLmks[0];
  const auto & mov_ref = movingLmks[0];
  for(size_t i = 1; i < numLmks; ++i)
  {
    double weight = ( has_lmk_weights) ? landmarkWgts[i] : 1.0;
    fix_dist_sum += weight*fix_ref.EuclideanDistanceTo(fixedLmks[i]);
    mov_dist_sum += weight*mov_ref.EuclideanDistanceTo(movingLmks[i]);
  }
  isotropic_scale_factor = fix_dist_sum/mov_dist_sum;
  return isotropic_scale_factor;
}

void
PreProcessLandmarkFiles(std::string              inputFixedLandmarkFilename,
                        std::string              inputMovingLandmarkFilename,
                        std::string              inputWeightFilename,
                        LandmarkPointContainer & fixedLmks,
                        LandmarkPointContainer & movingLmks,
                        LandmarkWeightType &     landmarkWgts
                        )
{
  fixedLmks.clear();
  movingLmks.clear();
  landmarkWgts.clear();


  /** read in *fcsv file */
  /** check four landmarks */
  std::cout << "Reading fixed landmarks set: " << inputFixedLandmarkFilename << std::endl;
  LandmarksMapType fixedLandmarks = ReadSlicer3toITKLmk(inputFixedLandmarkFilename);

  std::cout << "Reading moving landmarks set: " << inputMovingLandmarkFilename << std::endl;
  LandmarksMapType movingLandmarks = ReadSlicer3toITKLmk(inputMovingLandmarkFilename);

  /** Landmark Weights */
  LandmarksWeightMapType landmarkWeightMap;
  if (!inputWeightFilename.empty())
  {
    std::cout << "Reading landmarks weight file: " << inputWeightFilename << std::endl;
    landmarkWeightMap = ReadLandmarkWeights(inputWeightFilename);
    CheckLandmarks(fixedLandmarks, landmarkWeightMap);
    CheckLandmarks(movingLandmarks, landmarkWeightMap);
  }

  using LandmarkConstIterator = LandmarksMapType::const_iterator;
  for (LandmarkConstIterator fixedIt = fixedLandmarks.begin(); fixedIt != fixedLandmarks.end(); ++fixedIt)
  {
    LandmarkConstIterator movingIt = movingLandmarks.find(fixedIt->first);
    if (movingIt != movingLandmarks.end())
    {
      fixedLmks.push_back(fixedIt->second);
      movingLmks.push_back(movingIt->second);

      if (!landmarkWeightMap.empty())
      {
        if (landmarkWeightMap.find(fixedIt->first) != landmarkWeightMap.end())
        {
          landmarkWgts.push_back(landmarkWeightMap[fixedIt->first]);
        }
        else
        {
          std::cout << "Landmark for " << fixedIt->first << " does not exist. "
                    << "Set the weight to 0.5 " << std::endl;
          landmarkWgts.push_back(0.5F);
        }
      }
    }
  }
}


template <typename TTransformType>
typename TTransformType::Pointer
InitializeTransform(const LandmarkPointContainer & fixedLmks,
                    const LandmarkPointContainer & movingLmks,
                    LandmarkWeightType &           landmarkWgts,
                    ImageType::Pointer             referenceImage,
                    int                            bsplineNumberOfControlPoints)
{
  /** Landmark Initializaer */

  using LocalTransformType = TTransformType;
  typename LocalTransformType::Pointer transform = LocalTransformType::New();

  using LandmarkBasedInitializerType = itk::LandmarkBasedTransformInitializer<LocalTransformType, ImageType, ImageType>;

  typename LandmarkBasedInitializerType::Pointer landmarkBasedInitializer = LandmarkBasedInitializerType::New();


  /** set weights */
  if (!landmarkWgts.empty())
  {
    landmarkBasedInitializer->SetLandmarkWeight(landmarkWgts);
  }

  if (referenceImage.IsNotNull())
  {
    landmarkBasedInitializer->SetReferenceImage(referenceImage);
  }
  landmarkBasedInitializer->SetBSplineNumberOfControlPoints(bsplineNumberOfControlPoints);

  landmarkBasedInitializer->SetFixedLandmarks(fixedLmks);
  landmarkBasedInitializer->SetMovingLandmarks(movingLmks);
  landmarkBasedInitializer->SetTransform(transform);
  try
  {
    landmarkBasedInitializer->InitializeTransform();
  }
  catch (itk::ExceptionObject & err)
  {
    std::cout << "Exception Object caught: " << std::endl;
    std::cout << err << std::endl;
    throw;
  }
  return transform;
}

//////////////////// M A I N /////////////////////////////////////////////////
int
main(int argc, char * argv[])
{
  PARSE_ARGS;
  BRAINSRegisterAlternateIO();

  if (inputFixedLandmarkFilename.empty() || inputMovingLandmarkFilename.empty() || outputTransformFilename.empty())
  {
    std::cout << "Input Landmarks ( inputFixedLandmarkFilename ,"
              << "inputMovingLandmarkFilename ) and "
              << "outputTransformationFilename are necessary" << std::endl;
    exit(EXIT_FAILURE);
  }

  ImageType::Pointer referenceImage = nullptr;
  if (!inputReferenceImageFilename.empty())
  {
    using ReaderType = itk::ImageFileReader<ImageType>;
    ReaderType::Pointer reader = ReaderType::New();
    reader->SetFileName(inputReferenceImageFilename);
    try
    {
      reader->Update();
      referenceImage = reader->GetOutput();
    }
    catch (itk::ExceptionObject & err)
    {
      std::cout << "Exception Object caught: " << std::endl;
      std::cout << err << std::endl;
      throw;
    }
  }


  LandmarkPointContainer fixedLmks;
  LandmarkPointContainer movingLmks;
  LandmarkWeightType     landmarkWgts;

  PreProcessLandmarkFiles(
    inputFixedLandmarkFilename, inputMovingLandmarkFilename, inputWeightFilename, fixedLmks, movingLmks, landmarkWgts);

  // =================

  if (outputTransformType == "AffineTransform")
  {
    using AffineTransformType = itk::AffineTransform<ParameterValueType, Dimension>;

    auto transform = InitializeTransform<AffineTransformType>(
      fixedLmks, movingLmks, landmarkWgts, referenceImage, bsplineNumberOfControlPoints);
    std::cout << "Writing output transform file to disk: " << outputTransformFilename << " type:" << outputTransformType << std::endl;
    itk::WriteTransformToDisk<double>(transform, outputTransformFilename);
  }
  else if (outputTransformType == "BSplineTransform")
  {
    constexpr static unsigned int SplineOrder = 3;
    using BSplineTransformType = itk::BSplineTransform<ParameterValueType, Dimension, SplineOrder>;
    auto transform = InitializeTransform<BSplineTransformType>(
      fixedLmks, movingLmks, landmarkWgts, referenceImage, bsplineNumberOfControlPoints);
    std::cout << "Writing output transform file to disk: " << outputTransformFilename << std::endl;
    itk::WriteTransformToDisk<double>(transform, outputTransformFilename);
  }
  else if (outputTransformType == "VersorRigid3DTransform" || outputTransformType == "Similarity3DTransform")
  {
    using VersorRigid3DTransformType = itk::VersorRigid3DTransform<ParameterValueType>;
    auto transform = InitializeTransform<VersorRigid3DTransformType>(
      fixedLmks, movingLmks, landmarkWgts, referenceImage, bsplineNumberOfControlPoints);

    if (outputTransformType == "Similarity3DTransform")
    { // Now Add isotropic scaling
      using Similarity3DTransformType = itk::Similarity3DTransform<ParameterValueType>;
      Similarity3DTransformType::Pointer simTransform = Similarity3DTransformType::New();

      simTransform->SetMatrix(transform->GetMatrix());
      const double isotropic_scale_factor= ComputeIsotropicScaleFactor(fixedLmks, movingLmks, landmarkWgts);
      std::cout << "Using isotropic scale factor: " << isotropic_scale_factor << std::endl;
      //simTransform->SetScale(isotropic_scale_factor);
      std::cout << "Writing output transform file to disk: " << outputTransformFilename << " type:" << outputTransformType << std::endl;
      itk::WriteTransformToDisk<double>(simTransform, outputTransformFilename);
    }
    else
    {
      std::cout << "Writing output transform file to disk: " << outputTransformFilename << " type:" << outputTransformType << std::endl;
      itk::WriteTransformToDisk<double>(transform, outputTransformFilename);
    }
  }
  else
  {
    std::cerr << "Error: Invalid parameter for output transform type." << std::endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
