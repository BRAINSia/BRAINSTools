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

int main(int argc, char *argv[])
{
  // HACK:  NOTE THIS PROGRAM STILL USES ARGV ARGC, and ignores the PARSE_ARGS.
  // It needs to be fixed.
  PARSE_ARGS;
  BRAINSRegisterAlternateIO();
  // Image, filter typedef
  constexpr unsigned int LocalImageDimension = 3;

  typedef short PixelType;

  typedef itk::Image<PixelType, LocalImageDimension> ImageType;
  typedef ImageType::PointType                       ImagePointType;

  typedef itk::ImageFileReader<ImageType>         ReaderType;
  typedef itk::ImageFileWriter<ImageType>         WriterType;
  typedef itk::FindCenterOfBrainFilter<ImageType> FindCenterFilter;
  typedef itk::BRAINSHoughEyeDetector<
      ImageType, ImageType>                             HoughEyeDetectorType;

  // Read input image
  if( argc < 3 )
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
  catch( itk::ExceptionObject & err )
    {
    std::cerr << " Error while reading image file(s) with ITK:\n "
              << err << std::endl;
    }

  // Find center of head mass
  std::cout << "Finding center of head mass..." << std::endl;
  FindCenterFilter::Pointer findCenterFilter = FindCenterFilter::New();
  findCenterFilter->SetInput( reader->GetOutput() );
  findCenterFilter->SetAxis(2);
  findCenterFilter->SetOtsuPercentileThreshold(0.01);
  findCenterFilter->SetClosingSize(7);
  findCenterFilter->SetHeadSizeLimit(700);
  findCenterFilter->SetBackgroundValue(0);
  findCenterFilter->Update();
  ImagePointType centerOfHeadMass = findCenterFilter->GetCenterOfBrain();

  // Find eye centers with BRAINS Hough Eye Detector
  std::cout << "Finding eye centers..." << std::endl;
  HoughEyeDetectorType::Pointer houghEyeDetector = HoughEyeDetectorType::New();
  houghEyeDetector->SetInput( reader->GetOutput() );
  houghEyeDetector->SetHoughEyeDetectorMode(1); // For T1 images
  houghEyeDetector->SetCenterOfHeadMass(centerOfHeadMass);
  houghEyeDetector->SetResultsDir(debugDir);           // debug image write dir
  houghEyeDetector->SetWritedebuggingImagesLevel(2);   // write ROI and
                                                       // accumulator images
  houghEyeDetector->Update();

  ImagePointType leftEye   = houghEyeDetector->GetLE();
  ImagePointType rightEye  = houghEyeDetector->GetRE();
  ImagePointType alignedLE =
    houghEyeDetector->GetInvVersorTransform()->TransformPoint(leftEye);
  ImagePointType alignedRE =
    houghEyeDetector->GetInvVersorTransform()->TransformPoint(rightEye);

  std::cout << "Left eye = " << leftEye << std::endl;
  std::cout << "Right eye = " << rightEye << std::endl;
  std::cout << "Aligned left eye = " << alignedLE << std::endl;
  std::cout << "Alinged right eye = " << alignedRE << std::endl;

  // Write aligned image
  WriterType::Pointer writer = WriterType::New();
  writer->UseCompressionOn();
  writer->SetFileName(outputVolume);
  writer->SetInput( houghEyeDetector->GetOutput() );
  try
    {
    writer->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    std::cerr << " Error while writing image file(s) with ITK:\n "
              << err << std::endl;
    }
  return 0;
}
