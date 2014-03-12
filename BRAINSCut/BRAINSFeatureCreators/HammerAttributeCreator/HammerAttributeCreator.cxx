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
/******************************************************************************/
/*  Changed from Version 1.0                                                  */
/*  1. Subvolume matching is also added for subject driving voxels            */
/*      during correspondence detection.                                      */
/*  2. samll testing                                                          */
/*     (a) Binarize the tissue matching results  -> Good for Young Brain only */
/*     (b) The background effect in subvolume matching is also considered ->  */
/*                                                          NOT GOOD FOR ALL  */
/*  Search 'March ?? 2003' for these changes.                                 */
/******************************************************************************/
#include "itkImage.h"
#include "itkVector.h"
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#include "itkImageRegionIteratorWithIndex.h"

#define PARTIALVOLUME 1

#if PARTIALVOLUME > 0
#include "itkHammerTissueAttributeVectorFromPartialVolumeImageFilter.h"
#else
#include "itkHammerTissueAttributeVectorImageFilter.h"
#endif

#include "HammerAttributeCreatorCLP.h"
#include <BRAINSCommonLib.h>

typedef itk::Image<float, 3>               ImageType;
typedef itk::HammerTissueAttributeVector   AttributeVectorType;
typedef itk::Image<AttributeVectorType, 3> AttributeImageType;
#if PARTIALVOLUME  > 0
typedef itk::HammerTissueAttributeVectorFromPartialVolumeImageFilter<ImageType, AttributeImageType> AttributeFilterType;
#else
typedef itk::HammerTissueAttributeVectorImageFilter<ImageType, AttributeImageType> AttributeFilterType;
#endif

static void WriteAttributeComponent(const std::string& filename, AttributeImageType::Pointer avImg, const int idx)
{
  ImageType::Pointer cImg = ImageType::New();

  cImg->CopyInformation(avImg);
  cImg->SetRegions( cImg->GetLargestPossibleRegion() );
  cImg->Allocate();

  itk::ImageRegionIteratorWithIndex<ImageType> it( cImg, cImg->GetLargestPossibleRegion() );
  for( it.GoToBegin(); !it.IsAtEnd(); ++it )
    {
    ImageType::IndexType          imgIdx = it.GetIndex();
    AttributeImageType::PixelType p = avImg->GetPixel(imgIdx);
    it.Set(p[idx]);
    }

  itk::ImageFileWriter<ImageType>::Pointer w = itk::ImageFileWriter<ImageType>::New();
  w->UseCompressionOn();
  w->SetFileName(filename);
  w->SetInput(cImg);
  w->Update();
  return;
}

int main(int argc, char *argv[])
{
  PARSE_ARGS;
  BRAINSRegisterAlternateIO();
  /*** Load in fixed image and compute the attribute vectors ***/
  itk::ImageFileReader<ImageType>::Pointer gmReader
    = itk::ImageFileReader<ImageType>::New();
  itk::ImageFileReader<ImageType>::Pointer wmReader
    = itk::ImageFileReader<ImageType>::New();
  itk::ImageFileReader<ImageType>::Pointer csfReader
    = itk::ImageFileReader<ImageType>::New();
  gmReader->SetFileName(inputGMVolume);
  wmReader->SetFileName(inputWMVolume);
  csfReader->SetFileName(inputCSFVolume);
  try
    {
    gmReader->Update();
    wmReader->Update();
    csfReader->Update();
    }
  catch( itk::ExceptionObject *ex )
    {
    std::cerr << ex << std::endl;
    }
  ImageType::Pointer myImage = gmReader->GetOutput();

  // HACK: Need to replace some values.
  itk::ImageRegionIterator<ImageType> it( myImage, myImage->GetLargestPossibleRegion() );
  for( it.GoToBegin(); !it.IsAtEnd(); ++it )
    {
    const ImageType::PixelType tmp = it.Get();
    if( tmp == 3 )
      {
      it.Set(2);
      }
    if( tmp == 5 )
      {
      it.Set(0);
      }
    }

  std::cout << "Image read in\n";

  AttributeFilterType::Pointer modleAttributeFilter = AttributeFilterType::New();
  modleAttributeFilter->SetCSFVolume( wmReader->GetOutput() );
  modleAttributeFilter->SetGMVolume( gmReader->GetOutput() );
  modleAttributeFilter->SetWMVolume( wmReader->GetOutput() );

  // not thread-safe, yet!
  modleAttributeFilter->SetNumberOfThreads(1);
  modleAttributeFilter->SetStrength(Strength);
  modleAttributeFilter->SetScale(Scale);
  modleAttributeFilter->Update();

  AttributeImageType::Pointer fixedAVec = modleAttributeFilter->GetOutput();

  WriteAttributeComponent(outputVolumeBase + "_GMIEdgiInformation.nii.gz", fixedAVec, 0);
  WriteAttributeComponent(outputVolumeBase + "_GMIInputLabel.nii.gz", fixedAVec, 1);
  WriteAttributeComponent(outputVolumeBase + "_GMINonWM.nii.gz", fixedAVec, 2);
  WriteAttributeComponent(outputVolumeBase + "_GMICSF.nii.gz", fixedAVec, 3);
  WriteAttributeComponent(outputVolumeBase + "_GMIGM.nii.gz", fixedAVec, 4);

  return 0;
}
