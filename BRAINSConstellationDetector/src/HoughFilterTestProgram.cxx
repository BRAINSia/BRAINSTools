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
// Inspecting the performance of the HoughTransformRadialVotingImageFilter regarding the number of threads

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkHoughTransformRadialVotingImageFilter.h"

int
main(int argc, char * argv[])
{
  constexpr unsigned int LocalImageDimension = 3;

  using InputPixelType = short;
  using OutputPixelType = double;
  using InputImageType = itk::Image<InputPixelType, LocalImageDimension>;
  using OutputImageType = itk::Image<OutputPixelType, LocalImageDimension>;

  using HoughFilterType = itk::HoughTransformRadialVotingImageFilter<InputImageType, OutputImageType>;
  using SpheresListType = HoughFilterType::SpheresListType;
  using HoughFilterPointer = HoughFilterType::Pointer;
  HoughFilterPointer houghFilter = HoughFilterType::New();

  // Reader type
  using ReaderType = itk::ImageFileReader<InputImageType>;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName(argv[1]);
  // reader->Update();

  houghFilter->SetInput(reader->GetOutput());

  houghFilter->SetNumberOfSpheres(2);
  houghFilter->SetMinimumRadius(11.);
  houghFilter->SetMaximumRadius(13.);
  houghFilter->SetSigmaGradient(1.);
  houghFilter->SetVariance(1.);
  houghFilter->SetSphereRadiusRatio(1.);
  houghFilter->SetVotingRadiusRatio(.5);
  houghFilter->SetThreshold(10.);
  houghFilter->SetOutputThreshold(.8);
  houghFilter->SetGradientThreshold(0.);
  houghFilter->SetNbOfThreads(64);
  houghFilter->SetSamplingRatio(.2);
  houghFilter->SetHoughEyeDetectorMode(1);
  houghFilter->SetWritedebuggingAccumulatorImageLevel(0);

  try
  {
    houghFilter->Update();
  }
  catch (itk::ExceptionObject & excep)
  {
    std::cerr << "Failed houghFilter " << std::endl;
    std::cerr << excep << std::endl;
  }
  catch (...)
  {
    std::cout << "Failed on houghFilter exception occured" << std::endl;
  }

  /*
  this->m_AccumulatorImage = houghFilter->GetOutput();

  // Write debug accumulator image
  typename WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( this->m_ResultsDir + "/HoughEyeAccumulator.nii.gz" );
  writer->SetInput( this->m_AccumulatorImage );
  writer->SetUseCompression( true );

  */

  // writer type
  using WriterType = itk::ImageFileWriter<OutputImageType>;
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName(argv[2]);
  writer->SetInput(houghFilter->GetOutput());
  try
  {
    writer->Update();
  }
  catch (itk::ExceptionObject & excep)
  {
    std::cerr << "Cannot write the Accumulator image!" << std::endl;
    std::cerr << excep << std::endl;
  }

  const SpheresListType spheres = houghFilter->GetSpheres();
  if (spheres.size() < 2)
  {
    std::cerr << "Error: The number of detected spheres is less than 2!" << std::endl;
    // std::cerr << "The program will continue to run for generating some debug output for GUI corrector." << std::endl;
    // std::cerr << "Make sure you set the debug level > 4." << std::endl;
    // this->m_Failure = true;
    // return -1;
  }
  else
  {
    std::cout << "It works! The number of detected spheres is not less than 2!" << std::endl;
  }

  return EXIT_SUCCESS;
}
