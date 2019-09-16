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
/*=========================================================================

 Program:   GTRACT (Guided Tensor Restore Anatomical Connectivity Tractography)
 Module:    $RCSfile: $
 Language:  C++
 Date:      $Date: 2010/05/03 14:53:40 $
 Version:   $Revision: 1.9 $

   Copyright (c) University of Iowa Department of Radiology. All rights reserved.
   See GTRACT-Copyright.txt or http://mri.radiology.uiowa.edu/copyright/GTRACT-Copyright.txt
   for details.

      This software is distributed WITHOUT ANY WARRANTY; without even
      the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
      PURPOSE.  See the above copyright notices for more information.

=========================================================================*/

#include <iostream>
#include <fstream>

#include <itkImage.h>
#include <itkVectorIndexSelectionCastImageFilter.h>
#include <itkThresholdImageFilter.h>
#include <itkOrientImageFilter.h>

#include <itkCompositeTransform.h>
#include <itkTransform.h>

#include "itkGtractImageIO.h"
#include "BRAINSFitHelper.h"

#include "gtractCoRegAnatomyCLP.h"
#include "BRAINSThreadControl.h"
#include "DWIConvertLib.h"

int
main(int argc, char * argv[])
{
  PARSE_ARGS;
  BRAINSRegisterAlternateIO();
  const BRAINSUtils::StackPushITKDefaultNumberOfThreads TempDefaultNumberOfThreadsHolder(numberOfThreads);

  std::vector<int> GridSize;
  GridSize.push_back(gridSize[0]);
  GridSize.push_back(gridSize[1]);
  GridSize.push_back(gridSize[2]);

  bool debug = true;
  if (debug)
  {
    std::cout << "=====================================================" << std::endl;
    std::cout << "Input Image: " << inputVolume << std::endl;
    std::cout << "Output  Transform: " << outputTransformName << std::endl;
    std::cout << "Anatomical Image: " << inputAnatomicalVolume << std::endl;
    std::cout << "Iterations: " << numberOfIterations << std::endl;
    if (transformType == "Bspline")
    {
      std::cout << "Input Rigid Transform: " << inputRigidTransform << std::endl;
      // std::cout << "Grid Size: " << GridSize <<std::endl;
      // std::cout << "Border Size: " << borderSize <<std::endl;
      //    std::cout << "Corrections: " << numberOfCorrections <<std::endl;
      //    std::cout << "Evaluations: " << numberOfEvaluations <<std::endl;
      std::cout << "Histogram: " << numberOfHistogramBins << std::endl;
      std::cout << "Scale: " << spatialScale << std::endl;
      std::cout << "Convergence: " << convergence << std::endl;
      std::cout << "Gradient Tolerance: " << gradientTolerance << std::endl;
      std::cout << "Index: " << vectorIndex << std::endl;
    }
    else if (transformType == "Rigid")
    {
      std::cout << "Translation Scale: " << translationScale << std::endl;
      std::cout << "Maximum Step Length: " << maximumStepSize << std::endl;
      std::cout << "Minimum Step Length: " << minimumStepSize << std::endl;
      std::cout << "Relaxation Factor: " << relaxationFactor << std::endl;
      std::cout << "Samples: " << numberOfSamples << std::endl;
      std::cout << "Index: " << vectorIndex << std::endl;
    }
    //    std::cout << "Bound X: " << boundX <<std::endl;
    //    std::cout << "\tLower X Bound: " << xLowerBound <<std::endl;
    //    std::cout << "\tUpper X Bound: " << xUpperBound <<std::endl;
    //    std::cout << "Bound Y: " << boundY <<std::endl;
    //    std::cout << "\tLower Y Bound: " << yLowerBound <<std::endl;
    //    std::cout << "\tUpper Y Bound: " << yUpperBound <<std::endl;
    //    std::cout << "Bound Z: " << boundZ <<std::endl;
    //    std::cout << "\tLower Z Bound: " << zLowerBound <<std::endl;
    //    std::cout << "\tUpper Z Bound: " << zUpperBound <<std::endl;
    std::cout << "=====================================================" << std::endl;
  }

  bool violated = false;
  if (inputVolume.size() == 0)
  {
    violated = true;
    std::cout << "  --inputVolume Required! " << std::endl;
  }
  if (inputAnatomicalVolume.size() == 0)
  {
    violated = true;
    std::cout << "  --inputAnatomicalVolume Required! " << std::endl;
  }
  if (transformType == "Bspline")
  {
    if (inputRigidTransform.size() == 0)
    {
      violated = true;
      std::cout << "  --inputRigidTransform Required! " << std::endl;
    }
  }
  if (outputTransformName.size() == 0)
  {
    violated = true;
    std::cout << "  --outputTransform Required! " << std::endl;
  }
  if (violated)
  {
    return EXIT_FAILURE;
  }

  std::string convertedVolume;
  if (convertInputVolumeToNrrdOrNifti("Nrrd", inputVolume, convertedVolume))
  {
    inputVolume = convertedVolume;
  }
  else
  {
    std::cout << "Error: DWI Convert can not read inputVolume." << std::endl;
    return -1;
  }

  // using PixelType = signed short;
  using PixelType = float;
  using VectorImageType = itk::VectorImage<PixelType, 3>;

  using VectorImageReaderType = itk::ImageFileReader<VectorImageType, itk::DefaultConvertPixelTraits<PixelType>>;
  VectorImageReaderType::Pointer vectorImageReader = VectorImageReaderType::New();
  vectorImageReader->SetFileName(inputVolume);

  try
  {
    vectorImageReader->Update();
  }
  catch (itk::ExceptionObject & ex)
  {
    std::cout << ex << std::endl;
    throw;
  }

  using AnatomicalImageType = itk::Image<PixelType, 3>;
  using AnatomicalImageReaderType = itk::ImageFileReader<AnatomicalImageType>;
  AnatomicalImageReaderType::Pointer anatomicalReader = AnatomicalImageReaderType::New();
  anatomicalReader->SetFileName(inputAnatomicalVolume);

  try
  {
    anatomicalReader->Update();
  }
  catch (itk::ExceptionObject & ex)
  {
    std::cout << ex << std::endl;
    throw;
  }

  /* Extract the Vector Image Index for Registration */
  using VectorSelectFilterType = itk::VectorIndexSelectionCastImageFilter<VectorImageType, AnatomicalImageType>;
  using VectorSelectFilterPointer = VectorSelectFilterType::Pointer;

  VectorSelectFilterPointer selectIndexImageFilter = VectorSelectFilterType::New();
  selectIndexImageFilter->SetIndex(vectorIndex);
  selectIndexImageFilter->SetInput(vectorImageReader->GetOutput());
  try
  {
    selectIndexImageFilter->Update();
  }
  catch (itk::ExceptionObject & e)
  {
    std::cout << e << std::endl;
    throw;
  }

  std::string localInitializeTransformMode = "Off";
  if (((useCenterOfHeadAlign == true) + (useGeometryAlign == true) + (useMomentsAlign == true)) > 1)
  {
    std::cout << "ERROR:  Can only specify one of [useCenterOfHeadAlign | useGeometryAlign | useMomentsAlign ]"
              << std::endl;
  }
  if (useCenterOfHeadAlign == true)
  {
    localInitializeTransformMode = "useCenterOfHeadAlign";
  }
  if (useGeometryAlign == true)
  {
    localInitializeTransformMode = "useGeometryAlign";
  }
  if (useMomentsAlign == true)
  {
    localInitializeTransformMode = "useMomentsAlign";
  }

  using RegisterFilterType = itk::BRAINSFitHelper;
  RegisterFilterType::Pointer registerImageFilter = RegisterFilterType::New();

  if (transformType == "Rigid")
  {
    /* The Threshold Image Filter is used to produce the brain clipping mask. */
    using ThresholdFilterType = itk::ThresholdImageFilter<AnatomicalImageType>;
    constexpr PixelType          imageThresholdBelow = 100;
    ThresholdFilterType::Pointer brainOnlyFilter = ThresholdFilterType::New();
    brainOnlyFilter->SetInput(selectIndexImageFilter->GetOutput());
    brainOnlyFilter->ThresholdBelow(imageThresholdBelow);
    try
    {
      brainOnlyFilter->Update();
    }
    catch (itk::ExceptionObject & e)
    {
      std::cout << e << std::endl;
      throw;
    }
    registerImageFilter->SetMovingVolume(brainOnlyFilter->GetOutput());
  }
  if (transformType == "Bspline")
  {
    using OrientFilterType = itk::OrientImageFilter<AnatomicalImageType, AnatomicalImageType>;
    OrientFilterType::Pointer orientImageFilter = OrientFilterType::New();
    //  orientImageFilter->SetInput(brainOnlyFilter->GetOutput() );
    orientImageFilter->SetInput(selectIndexImageFilter->GetOutput());
    orientImageFilter->SetDesiredCoordinateDirection(anatomicalReader->GetOutput()->GetDirection());
    orientImageFilter->UseImageDirectionOn();
    try
    {
      orientImageFilter->Update();
    }
    catch (itk::ExceptionObject & e)
    {
      std::cout << e << std::endl;
      throw;
    }
    registerImageFilter->SetMovingVolume(orientImageFilter->GetOutput());
  }

  std::vector<std::string> transformTypes;
  std::vector<int>         iterations;
  iterations.push_back(numberOfIterations);

  using TransformType = itk::Transform<double, 3, 3>;
  using CompositeTransformType = itk::CompositeTransform<double, 3>;

  if (transformType == "Bspline")
  {
    transformTypes.push_back("BSpline");

    registerImageFilter->SetSamplingPercentage(1.0 / spatialScale);
    registerImageFilter->SetNumberOfHistogramBins(numberOfHistogramBins);
    registerImageFilter->SetSplineGridSize(gridSize);
    registerImageFilter->SetCostFunctionConvergenceFactor(convergence);
    registerImageFilter->SetProjectedGradientTolerance(gradientTolerance);
    registerImageFilter->SetMaxBSplineDisplacement(maxBSplineDisplacement);
    registerImageFilter->SetInitializeTransformMode(localInitializeTransformMode);
    if (inputRigidTransform.size() > 0)
    {
      TransformType::Pointer          inputTransform = itk::ReadTransformFromDisk(inputRigidTransform);
      CompositeTransformType::Pointer inputCompositeTransform =
        dynamic_cast<CompositeTransformType *>(inputTransform.GetPointer());
      if (inputCompositeTransform.IsNull())
      {
        inputCompositeTransform = CompositeTransformType::New();
        inputCompositeTransform->AddTransform(inputTransform);
      }
      registerImageFilter->SetCurrentGenericTransform(inputCompositeTransform);
    }
  }

  if (transformType == "Rigid")
  {
    transformTypes.push_back("ScaleVersor3D");
    std::vector<double> minStepLength;
    minStepLength.push_back((double)minimumStepSize);
    registerImageFilter->SetTranslationScale(translationScale);
    registerImageFilter->SetMaximumStepLength(maximumStepSize);
    registerImageFilter->SetMinimumStepLength(minStepLength);
    registerImageFilter->SetRelaxationFactor(relaxationFactor);
    if (numberOfSamples > 0)
    {
      const unsigned long numberOfAllSamples = anatomicalReader->GetOutput()->GetBufferedRegion().GetNumberOfPixels();
      samplingPercentage = static_cast<double>(numberOfSamples) / numberOfAllSamples;
      std::cout << "WARNING --numberOfSamples is deprecated, please use --samplingPercentage instead " << std::endl;
      std::cout << "WARNING: Replacing command line --samplingPercentage " << samplingPercentage << std::endl;
    }
    registerImageFilter->SetSamplingPercentage(samplingPercentage);
    registerImageFilter->SetInitializeTransformMode(localInitializeTransformMode);
  }
  registerImageFilter->SetFixedVolume(anatomicalReader->GetOutput());
  registerImageFilter->SetTransformType(transformTypes);
  registerImageFilter->SetNumberOfIterations(iterations);
  try
  {
    registerImageFilter->Update();
  }
  catch (itk::ExceptionObject & ex)
  {
    std::cout << ex << std::endl;
    throw;
  }
  using GenericTransformType = itk::Transform<double, 3, 3>;
  GenericTransformType::Pointer outputTransform = registerImageFilter->GetCurrentGenericTransform()->GetNthTransform(0);
  itk::WriteTransformToDisk<double>(outputTransform, outputTransformName);
  return EXIT_SUCCESS;
}
