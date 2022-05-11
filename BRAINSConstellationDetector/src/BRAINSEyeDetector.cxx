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
 * at Psychiatry Imaging Lab,
 * University of Iowa Health Care 2010
 */

#include "itkFindCenterOfBrainFilter.h"
#include "BRAINSHoughEyeDetector.h"
#include "BRAINSEyeDetectorCLP.h"

#include "itkCommand.h"
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#include <iostream>
#include <BRAINSCommonLib.h>

int
main(int argc, char * argv[])
{
  // HACK:  NOTE THIS PROGRAM STILL USES ARGV ARGC, and ignores the PARSE_ARGS.
  // It needs to be fixed.
  PARSE_ARGS;
  BRAINSRegisterAlternateIO();
  // Image, filter typedef
  constexpr unsigned int LocalImageDimension = 3;

  using PixelType = short;

  using ImageType = itk::Image<PixelType, LocalImageDimension>;
  using ImagePointType = ImageType::PointType;

  using ReaderType = itk::ImageFileReader<ImageType>;
  using WriterType = itk::ImageFileWriter<ImageType>;
  using FindCenterFilter = itk::FindCenterOfBrainFilter<ImageType>;
  using HoughEyeDetectorType = itk::BRAINSHoughEyeDetector<ImageType, ImageType>;

  // Read input image
  if (argc < 3)
  {
    std::cerr << "Please specify the input filename." << std::endl;
    return EXIT_FAILURE;
  }

  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName(inputVolume);
  try
  {
    reader->Update();
  }
  catch (const itk::ExceptionObject & err)
  {
    std::cerr << " Error while reading image file(s) with ITK:\n " << err << std::endl;
  }

  // Find center of head mass
  std::cout << "Finding center of head mass..." << std::endl;
  FindCenterFilter::Pointer findCenterFilter = FindCenterFilter::New();
  findCenterFilter->SetInput(reader->GetOutput());
  findCenterFilter->SetAxis(2);
  findCenterFilter->SetOtsuPercentileThreshold(0.01);
  findCenterFilter->SetClosingSize(7);
  findCenterFilter->SetHeadSizeLimit(700);
  findCenterFilter->SetBackgroundValue(0);
  findCenterFilter->Update();
  ImagePointType orig_lmk_CenterOfHeadMass = findCenterFilter->GetCenterOfBrain();

  // Find eye centers with BRAINS Hough Eye Detector
  std::cout << "Finding eye centers..." << std::endl;
  HoughEyeDetectorType::Pointer houghEyeDetector = HoughEyeDetectorType::New();
  houghEyeDetector->SetInput(reader->GetOutput());
  houghEyeDetector->SetHoughEyeDetectorMode(1); // For T1 images
  houghEyeDetector->Setorig_lmk_CenterOfHeadMass(orig_lmk_CenterOfHeadMass);
  houghEyeDetector->SetResultsDir(debugDir);         // debug image write dir
  houghEyeDetector->SetWritedebuggingImagesLevel(2); // write ROI and
                                                     // accumulator images
  houghEyeDetector->Update();

  ImagePointType leftEye = houghEyeDetector->Getorig_lmk_LE();
  ImagePointType rightEye = houghEyeDetector->Getorig_lmk_RE();

  itk::VersorRigid3DTransform<double>::Pointer orig2eyeFixed_img_tfm =
    houghEyeDetector->GetModifiableorig2eyeFixedTransform();
  itk::VersorRigid3DTransform<double>::Pointer orig2eyeFixed_lmk_tfm = itk::VersorRigid3DTransform<double>::New();
  orig2eyeFixed_img_tfm->GetInverse(orig2eyeFixed_lmk_tfm);

  ImagePointType alignedLE = orig2eyeFixed_lmk_tfm->TransformPoint(leftEye);
  ImagePointType alignedRE = orig2eyeFixed_lmk_tfm->TransformPoint(rightEye);

  std::cout << "Left eye = " << leftEye << std::endl;
  std::cout << "Right eye = " << rightEye << std::endl;
  std::cout << "Aligned left eye = " << alignedLE << std::endl;
  std::cout << "Aligned right eye = " << alignedRE << std::endl;

  // Write aligned image
  WriterType::Pointer writer = WriterType::New();
  writer->UseCompressionOn();
  writer->SetFileName(outputVolume);
  writer->SetInput(houghEyeDetector->GetOutput());
  try
  {
    writer->Update();
  }
  catch (const itk::ExceptionObject & err)
  {
    std::cerr << " Error while writing image file(s) with ITK:\n " << err << std::endl;
  }
  return 0;
}
