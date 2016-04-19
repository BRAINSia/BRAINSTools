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
int main(int argc, char *argv[])
{
  PARSE_ARGS;
  BRAINSRegisterAlternateIO();

  typedef itk::Image<float, 3>         ProbabilityImageType;
  typedef double                       FloatingPointPrecision;

  typedef std::vector<ProbabilityImageType::Pointer> ProbabilityImageVectorType;
  typedef std::vector<bool>                          BoolVectorType;
  typedef vnl_vector<unsigned int>                   UnsignedIntVectorType;

  if( inputProbabilityVolume.size() < 1 )
    {
    std::cerr << "Missing probability volume list" << std::endl;
    return 1;
    }

  ProbabilityImageVectorType Posteriors;
  for( unsigned i = 0; i < inputProbabilityVolume.size(); ++i )
    {
    typedef itk::ImageFileReader<ProbabilityImageType> ImageReaderType;
    ImageReaderType::Pointer reader = ImageReaderType::New();
    reader->SetFileName(inputProbabilityVolume[i]);
    ProbabilityImageType::Pointer current;
    try
      {
      reader->Update();
      current = reader->GetOutput();
      }
    catch( itk::ExceptionObject & err )
      {
      std::cerr << err << " " << __FILE__ << " " << __LINE__ << std::endl;
      return 1;
      }
    Posteriors.push_back(current);
    }

  if( priorLabelCodes.size() < 1 )
    {
    std::cerr << "Missing prior label codes" << std::endl;
    }

  UnsignedIntVectorType priorLabels;
  priorLabels.set_size(priorLabelCodes.size() );
  for( unsigned i = 0; i < priorLabelCodes.size(); ++i )
    {
    priorLabels[i] = priorLabelCodes[i];
    }

  if( foregroundPriors.size() < 1 )
    {
    std::cerr << "Missing prior label codes" << std::endl;
    }

  BoolVectorType priorIsForeground;
  for( unsigned i = 0; i < foregroundPriors.size(); ++i )
    {
    priorIsForeground.push_back(foregroundPriors[i]);
    }

  ByteImageType::Pointer nonAirVolume;
  if( nonAirRegionMask == "" )
    {
    nonAirVolume = ByteImageType::New();
    ByteImageType::RegionType region =
      Posteriors[0]->GetLargestPossibleRegion();
    nonAirVolume->CopyInformation(Posteriors[0]);
    nonAirVolume->SetRegions(region);
    nonAirVolume->Allocate();
    ByteImageType::PixelType one = 1;
    nonAirVolume->FillBuffer(one);
    }
  else
    {
    typedef itk::ImageFileReader<ByteImageType> ImageReaderType;
    ImageReaderType::Pointer reader = ImageReaderType::New();
    reader->SetFileName(nonAirRegionMask);
    try
      {
      reader->Update();
      nonAirVolume = reader->GetOutput();
      }
    catch( itk::ExceptionObject & err )
      {
      std::cerr << err << " " << __FILE__ << " " << __LINE__ << std::endl;
      return 1;
      }
    }

  ByteImageType::Pointer dirtyLabels;
  ByteImageType::Pointer cleanLabels;

  try
    {
    ComputeLabels<ProbabilityImageType,
                  ByteImageType,
                  FloatingPointPrecision>
      (Posteriors,
      priorIsForeground,
      priorLabels,
      nonAirVolume,
      dirtyLabels,
      cleanLabels,
      inclusionThreshold, 0);
    }
  catch( itk::ExceptionObject & err )
    {
    std::cerr << err << " " << __FILE__ << " " << __LINE__ << std::endl;
    return 1;
    }

  if( dirtyLabelVolume != "" )
    {
    try
      {
      itkUtil::WriteImage<ByteImageType>(dirtyLabels, dirtyLabelVolume);
      }
    catch( itk::ExceptionObject & err )
      {
      std::cerr << err << " " << __FILE__ << " " << __LINE__ << std::endl;
      return 1;
      }
    }
  if( cleanLabelVolume != "" )
    {
    try
      {
      itkUtil::WriteImage<ByteImageType>(cleanLabels, cleanLabelVolume);
      }
    catch( itk::ExceptionObject & err )
      {
      std::cerr << err << " " << __FILE__ << " " << __LINE__ << std::endl;
      return 1;
      }
    }
  return 0;
}
