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
#include "BRAINSCreateLabelMapFromProbabilityMapsCLP.h"
#include "BRAINSComputeLabels.h"
#include "BRAINSCommonLib.h"
#include "itkIO.h"
#include "itkImage.h"
#include "itkImageFileReader.h"
int
main(int argc, char * argv[])
{
  PARSE_ARGS;
  BRAINSRegisterAlternateIO();

  using ProbabilityImageType = itk::Image<float, 3>;
  using FloatingPointPrecision = double;

  using ProbabilityImageVectorType = std::vector<ProbabilityImageType::Pointer>;
  using BoolVectorType = std::vector<bool>;
  using UnsignedIntVectorType = vnl_vector<unsigned int>;

  if (inputProbabilityVolume.empty())
  {
    std::cerr << "Missing probability volume list" << std::endl;
    return 1;
  }

  ProbabilityImageVectorType Posteriors;
  for (auto & i : inputProbabilityVolume)
  {
    using ImageReaderType = itk::ImageFileReader<ProbabilityImageType>;
    ImageReaderType::Pointer reader = ImageReaderType::New();
    reader->SetFileName(i);
    ProbabilityImageType::Pointer current;
    try
    {
      reader->Update();
      current = reader->GetOutput();
    }
    catch (const itk::ExceptionObject & err)
    {
      std::cerr << err << " " << __FILE__ << " " << __LINE__ << std::endl;
      return 1;
    }
    Posteriors.push_back(current);
  }

  if (priorLabelCodes.empty())
  {
    std::cerr << "Missing prior label codes" << std::endl;
  }

  UnsignedIntVectorType priorLabels;
  priorLabels.set_size(priorLabelCodes.size());
  for (unsigned i = 0; i < priorLabelCodes.size(); ++i)
  {
    priorLabels[i] = priorLabelCodes[i];
  }

  if (foregroundPriors.empty())
  {
    std::cerr << "Missing prior label codes" << std::endl;
  }

  BoolVectorType priorIsForeground;
  for (int foregroundPrior : foregroundPriors)
  {
    priorIsForeground.push_back(foregroundPrior);
  }

  ByteImageType::Pointer nonAirVolume;
  if (nonAirRegionMask.empty())
  {
    nonAirVolume = ByteImageType::New();
    ByteImageType::RegionType region = Posteriors[0]->GetLargestPossibleRegion();
    nonAirVolume->CopyInformation(Posteriors[0]);
    nonAirVolume->SetRegions(region);
    nonAirVolume->Allocate();
    ByteImageType::PixelType one = 1;
    nonAirVolume->FillBuffer(one);
  }
  else
  {
    using ImageReaderType = itk::ImageFileReader<ByteImageType>;
    ImageReaderType::Pointer reader = ImageReaderType::New();
    reader->SetFileName(nonAirRegionMask);
    try
    {
      reader->Update();
      nonAirVolume = reader->GetOutput();
    }
    catch (const itk::ExceptionObject & err)
    {
      std::cerr << err << " " << __FILE__ << " " << __LINE__ << std::endl;
      return 1;
    }
  }

  ByteImageType::Pointer dirtyLabels;
  ByteImageType::Pointer cleanLabels;

  try
  {
    ComputeLabels<ProbabilityImageType, ByteImageType, FloatingPointPrecision>(
      Posteriors, priorIsForeground, priorLabels, nonAirVolume, dirtyLabels, cleanLabels, inclusionThreshold, 0);
  }
  catch (const itk::ExceptionObject & err)
  {
    std::cerr << err << " " << __FILE__ << " " << __LINE__ << std::endl;
    return 1;
  }

  if (!dirtyLabelVolume.empty())
  {
    try
    {
      itkUtil::WriteImage<ByteImageType>(dirtyLabels, dirtyLabelVolume);
    }
    catch (const itk::ExceptionObject & err)
    {
      std::cerr << err << " " << __FILE__ << " " << __LINE__ << std::endl;
      return 1;
    }
  }
  if (!cleanLabelVolume.empty())
  {
    try
    {
      itkUtil::WriteImage<ByteImageType>(cleanLabels, cleanLabelVolume);
    }
    catch (const itk::ExceptionObject & err)
    {
      std::cerr << err << " " << __FILE__ << " " << __LINE__ << std::endl;
      return 1;
    }
  }
  return 0;
}
