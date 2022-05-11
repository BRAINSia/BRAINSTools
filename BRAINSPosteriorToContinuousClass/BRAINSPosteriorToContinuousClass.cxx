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
#include <iostream>
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionConstIterator.h"
#include "itkImageRegionIterator.h"
#include "BRAINSPosteriorToContinuousClassCLP.h"
#include "BRAINSCommonLib.h"

int
main(int argc, char ** argv)
{
  PARSE_ARGS;
  BRAINSRegisterAlternateIO();

  bool violated = false;
  if (inputWhiteVolume.empty())
  {
    violated = true;
    std::cout << "  --inputWhiteVolume Required! " << std::endl;
  }
  if (inputBasalGmVolume.empty())
  {
    violated = true;
    std::cout << "  --inputBasalGmVolume Required! " << std::endl;
  }
  if (inputSurfaceGmVolume.empty())
  {
    violated = true;
    std::cout << "  --inputSurfaceGmVolume Required! " << std::endl;
  }
  if (inputCsfVolume.empty())
  {
    violated = true;
    std::cout << "  --inputCsfVolume Required! " << std::endl;
  }
  if (inputVbVolume.empty())
  {
    violated = true;
    std::cout << "  --inputVbVolume Required! " << std::endl;
  }
  if (inputCrblGmVolume.empty())
  {
    violated = true;
    std::cout << "  --inputCrblGmVolume Required! " << std::endl;
  }
  if (inputCrblWmVolume.empty())
  {
    violated = true;
    std::cout << "  --inputCrblWmVolume Required! " << std::endl;
  }
  if (outputVolume.empty())
  {
    violated = true;
    std::cout << "  --outputVolume Required! " << std::endl;
  }
  if (violated)
  {
    exit(1);
  }

  using PixelType = float;
  using OutputPixelType = unsigned char;
  constexpr unsigned int Dimension = 3;

  using ImageType = itk::Image<PixelType, Dimension>;
  using OutputImageType = itk::Image<OutputPixelType, Dimension>;
  using ReaderType = itk::ImageFileReader<ImageType>;
  using WriterType = itk::ImageFileWriter<OutputImageType>;
  using ImageRegionConstIteratorType = itk::ImageRegionConstIterator<ImageType>;
  using ImageRegionIteratorType = itk::ImageRegionIterator<OutputImageType>;

  ReaderType::Pointer wmReader = ReaderType::New();
  ReaderType::Pointer basalGmReader = ReaderType::New();
  ReaderType::Pointer surfaceGmReader = ReaderType::New();
  ReaderType::Pointer csfReader = ReaderType::New();
  ReaderType::Pointer vbReader = ReaderType::New();
  ReaderType::Pointer crblGmReader = ReaderType::New();
  ReaderType::Pointer crblWmReader = ReaderType::New();
  WriterType::Pointer outputWriter = WriterType::New();

  ImageType::Pointer       wmVolume;
  ImageType::Pointer       bgmVolume;
  ImageType::Pointer       sgmVolume;
  ImageType::Pointer       crblGmVolume;
  ImageType::Pointer       crblWmVolume;
  ImageType::Pointer       csfVolume;
  ImageType::Pointer       vbVolume;
  OutputImageType::Pointer classVolume = OutputImageType::New();

  try
  {
    wmReader->SetFileName(inputWhiteVolume);
    wmReader->Update();
    wmVolume = wmReader->GetOutput();

    basalGmReader->SetFileName(inputBasalGmVolume);
    basalGmReader->Update();
    bgmVolume = basalGmReader->GetOutput();

    surfaceGmReader->SetFileName(inputSurfaceGmVolume);
    surfaceGmReader->Update();
    sgmVolume = surfaceGmReader->GetOutput();

    csfReader->SetFileName(inputCsfVolume);
    csfReader->Update();
    csfVolume = csfReader->GetOutput();

    vbReader->SetFileName(inputVbVolume);
    vbReader->Update();
    vbVolume = vbReader->GetOutput();

    crblGmReader->SetFileName(inputCrblGmVolume);
    crblGmReader->Update();
    crblGmVolume = crblGmReader->GetOutput();

    crblWmReader->SetFileName(inputCrblWmVolume);
    crblWmReader->Update();
    crblWmVolume = crblWmReader->GetOutput();
  }
  catch (const itk::ExceptionObject & exe)
  {
    std::cout << exe << std::endl;
    exit(1);
  }

  /* Allocate Class Voilume */
  classVolume->SetDirection(wmVolume->GetDirection());
  classVolume->SetRegions(wmVolume->GetLargestPossibleRegion());
  classVolume->SetOrigin(wmVolume->GetOrigin());
  classVolume->SetSpacing(wmVolume->GetSpacing());
  classVolume->Allocate();
  classVolume->FillBuffer(0);

  ImageRegionConstIteratorType wmItr(wmVolume, wmVolume->GetRequestedRegion());
  ImageRegionConstIteratorType bgmItr(bgmVolume, bgmVolume->GetRequestedRegion());
  ImageRegionConstIteratorType sgmItr(sgmVolume, sgmVolume->GetRequestedRegion());
  ImageRegionConstIteratorType csfItr(csfVolume, csfVolume->GetRequestedRegion());
  ImageRegionConstIteratorType vbItr(vbVolume, vbVolume->GetRequestedRegion());
  ImageRegionConstIteratorType crblGmItr(crblGmVolume, crblGmVolume->GetRequestedRegion());
  ImageRegionConstIteratorType crblWmItr(crblWmVolume, crblWmVolume->GetRequestedRegion());
  ImageRegionIteratorType      outItr(classVolume, classVolume->GetRequestedRegion());
  wmItr.GoToBegin();
  bgmItr.GoToBegin();
  sgmItr.GoToBegin();
  csfItr.GoToBegin();
  vbItr.GoToBegin();
  crblGmItr.GoToBegin();
  crblWmItr.GoToBegin();

  float minThreshold = 0.2;
  for (outItr.GoToBegin(); !outItr.IsAtEnd(); ++outItr)
  {
    float         maxGm = fmax(bgmItr.Value(), fmax(sgmItr.Value(), crblGmItr.Value()));
    float         maxWm = fmax(wmItr.Value(), crblWmItr.Value());
    float         maxCsf = csfItr.Value();
    float         maxVb = vbItr.Value();
    float         total;
    unsigned char voxelValue;

    if ((maxWm > maxGm) && (maxWm > maxCsf) && (maxWm > maxVb) && (maxWm > minThreshold))
    {
      total = maxWm + maxGm;
      voxelValue = itk::Math::rnd(130.0 + 120.0 * maxWm / total);
    }
    else if ((maxGm >= maxWm) && (maxGm >= maxCsf) && (maxGm > maxVb) && (maxGm > minThreshold))
    {
      if (maxWm >= maxCsf)
      {
        total = maxWm + maxGm;
        voxelValue = itk::Math::rnd(130.0 + 120.0 * maxWm / total);
      }
      else
      {
        total = maxCsf + maxGm;
        voxelValue = itk::Math::rnd(130.0 - 120.0 * maxCsf / total);
      }
    }
    else if ((maxCsf > maxGm) && (maxCsf >= maxWm) && (maxCsf > maxVb) && (maxCsf > minThreshold))
    {
      total = maxCsf + maxGm;
      voxelValue = itk::Math::rnd(10.0 + 120.0 * maxGm / total);
    }
    else if (maxVb > minThreshold)
    {
      voxelValue = 1;
    }
    else
    {
      voxelValue = 0;
    }
    outItr.Set(voxelValue);
    ++wmItr;
    ++bgmItr;
    ++sgmItr;
    ++csfItr;
    ++vbItr;
    ++crblGmItr;
    ++crblWmItr;
  }

  outputWriter->SetInput(classVolume);
  outputWriter->SetFileName(outputVolume);
  outputWriter->Update();

  return EXIT_SUCCESS;
}
