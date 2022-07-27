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
 Date:      $Date: 2006/03/29 14:53:40 $
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
#include <itkMedianImageFilter.h>
#include <itkResampleImageFilter.h>
#include <itkLinearInterpolateImageFunction.h>
#include <itkImageFileWriter.h>
#include <itkImageFileReader.h>
#include <itkComposeImageFilter.h>

#include <itkExtractImageFilter.h>
#include <itkImageMaskSpatialObject.h>

#include <itkDiffusionTensor3DReconstructionWithMaskImageFilter.h>
#include "itkGtractImageIO.h"
#include "itkGtractParameterIO.h"
#include "itkComputeDiffusionTensorImageFilter.h"
#include "itkResampleImageFilter.h"
#include "itkIdentityTransform.h"
#include "itkMetaDataObject.h"
#include "itkIOCommon.h"
#include "gtractTensorCLP.h"
#include "BRAINSThreadControl.h"

#include "GenericTransformImage.h"
#include "itkBRAINSROIAutoImageFilter.h"
#include "itkIO.h"

int
main(int argc, char * argv[])
{
  PARSE_ARGS;
  BRAINSRegisterAlternateIO();
  const BRAINSUtils::StackPushITKDefaultNumberOfThreads TempDefaultNumberOfThreadsHolder(numberOfThreads);

  using PixelType = signed short;
  using TensorPixelType = double;
  using VectorImageType = itk::VectorImage<PixelType, 3>;
  using IndexImageType = itk::Image<PixelType, 3>;
  using MaskImageType = itk::Image<unsigned char, 3>;
  IndexImageType::SizeType MedianFilterSize;
  MedianFilterSize[0] = medianFilterSize[0];
  MedianFilterSize[1] = medianFilterSize[1];
  MedianFilterSize[2] = medianFilterSize[2];

  bool debug = true;
  // applyMeasurementFrame = true;
  if (debug)
  {
    std::cout << "=====================================================" << std::endl;
    std::cout << "Input Image: " << inputVolume << std::endl;
    std::cout << "Output Image: " << outputVolume << std::endl;
    std::cout << "Resample Isotropic: " << resampleIsotropic << std::endl;
    std::cout << "Voxel Size: " << voxelSize << std::endl;
    std::cout << "Median Filter Size: " << MedianFilterSize << std::endl;
    std::cout << "Threshold: " << backgroundSuppressingThreshold << std::endl;
    std::cout << "B0 Index: " << b0Index << std::endl;
    std::cout << "Apply Measurement Frame: " << applyMeasurementFrame << std::endl;
    std::cout << "=====================================================" << std::endl;
  }

  bool violated = false;
  if (inputVolume.empty())
  {
    violated = true;
    std::cout << "  --inputVolume Required! " << std::endl;
  }
  if (outputVolume.empty())
  {
    violated = true;
    std::cout << "  --outputVolume Required! " << std::endl;
  }
  if (violated)
  {
    return EXIT_FAILURE;
  }

  using VectorImageReaderType = itk::ImageFileReader<VectorImageType, itk::DefaultConvertPixelTraits<PixelType>>;
  VectorImageReaderType::Pointer vectorImageReader = VectorImageReaderType::New();
  vectorImageReader->SetFileName(inputVolume);

  try
  {
    vectorImageReader->Update();
  }
  catch (const itk::ExceptionObject & ex)
  {
    std::cout << ex << std::endl;
    throw;
  }
  // figure out mask processing
  MaskImageType::ConstPointer maskImage; // will stay NULL if no mask is used.

  if (maskProcessingMode == "ROIAUTO")
  {
    using LocalIndexImageType = itk::Image<PixelType, 3>;
    using VectorSelectFilterType = itk::VectorIndexSelectionCastImageFilter<VectorImageType, LocalIndexImageType>;
    VectorSelectFilterType::Pointer SelectIndexImageFilter = VectorSelectFilterType::New();
    SelectIndexImageFilter->SetIndex(b0Index);
    SelectIndexImageFilter->SetInput(vectorImageReader->GetOutput());
    try
    {
      SelectIndexImageFilter->Update();
    }
    catch (const itk::ExceptionObject & e)
    {
      std::cerr << e << std::endl;
      return EXIT_FAILURE;
    };
    LocalIndexImageType::Pointer b0Image(SelectIndexImageFilter->GetOutput());
    using ROIAutoType = itk::BRAINSROIAutoImageFilter<LocalIndexImageType, MaskImageType>;
    ROIAutoType::Pointer ROIFilter = ROIAutoType::New();
    ROIFilter->SetInput(b0Image);
    ROIFilter->SetClosingSize(9.0); // default TODO Make parameter
    ROIFilter->SetDilateSize(0.0);
    ROIFilter->Update();
    maskImage = ROIFilter->GetBinaryImageROI();
  }
  else if (maskProcessingMode == "ROI")
  {
    if (maskVolume.empty())
    {
      std::cerr << "Error: missing mask Volume needed for ROI mask Processing" << std::endl;
      return EXIT_FAILURE;
    }
    maskImage = itkUtil::ReadImage<MaskImageType>(maskVolume);
    if (maskImage.IsNull())
    {
      std::cerr << "Error: can't read mask volume " << maskVolume << std::endl;
      return EXIT_FAILURE;
    }
  }
  /* Extract Diffusion Information from the Header */
  std::string BValue_str;
  std::string BValue_keyStr("DWMRI_b-value");
  itk::ExposeMetaData<std::string>(
    vectorImageReader->GetOutput()->GetMetaDataDictionary(), BValue_keyStr.c_str(), BValue_str);
  double BValue = std::stod(BValue_str.c_str());
  std::cout << "The BValue was found to be " << BValue_str << std::endl;

  std::vector<std::vector<double>> msrFrame;
  itk::ExposeMetaData<std::vector<std::vector<double>>>(
    vectorImageReader->GetOutput()->GetMetaDataDictionary(), "NRRD_measurement frame", msrFrame);
  TMatrix measurementFrame(3, 3);
  for (int i = 0; i < 3; i++)
  {
    for (int j = 0; j < 3; j++)
    {
      measurementFrame[i][j] = msrFrame[i][j];
    }
  }

  /* Process Invidual B-value Images and Reassemble the Vector Image */
  using TensorFilterType =
    itk::DiffusionTensor3DReconstructionWithMaskImageFilter<PixelType, PixelType, TensorPixelType>;
  using DirectionContainerType = TensorFilterType::GradientDirectionContainerType;
  DirectionContainerType::Pointer gradientDirectionContainer = DirectionContainerType::New();

  using VectorImageFilterType = itk::ComposeImageFilter<IndexImageType>;
  VectorImageFilterType::Pointer indexImageToVectorImageFilter = VectorImageFilterType::New();
  int                            vectorIndex = 0;
  for (unsigned int i = 0; i < vectorImageReader->GetOutput()->GetVectorLength(); i++)
  {
    using VectorSelectFilterType = itk::VectorIndexSelectionCastImageFilter<VectorImageType, IndexImageType>;
    using VectorSelectFilterPointer = VectorSelectFilterType::Pointer;

    VectorSelectFilterPointer selectIndexImageFilter = VectorSelectFilterType::New();
    selectIndexImageFilter->SetIndex(i);
    selectIndexImageFilter->SetInput(vectorImageReader->GetOutput());
    try
    {
      selectIndexImageFilter->Update();
    }
    catch (const itk::ExceptionObject & e)
    {
      std::cout << e << std::endl;
    }

    /* Median Filter */
    IndexImageType::Pointer baseImage;
    if (MedianFilterSize[0] > 0 || MedianFilterSize[1] > 0 || MedianFilterSize[2] > 0)
    {
      using MedianFilterType = itk::MedianImageFilter<IndexImageType, IndexImageType>;
      MedianFilterType::Pointer filter = MedianFilterType::New();
      filter->SetInput(selectIndexImageFilter->GetOutput());
      filter->SetRadius(MedianFilterSize);
      filter->Update();
      baseImage = filter->GetOutput();
    }
    else
    {
      baseImage = selectIndexImageFilter->GetOutput();
    }

    /* Resample To Isotropic Images */
    IndexImageType::Pointer bvalueImage;
    if (resampleIsotropic)
    {
      using ResampleFilterType = itk::ResampleImageFilter<IndexImageType, IndexImageType>;
      ResampleFilterType::Pointer resampler = ResampleFilterType::New();
      resampler->SetInput(baseImage);

      using InterpolatorType = itk::LinearInterpolateImageFunction<IndexImageType, double>;
      InterpolatorType::Pointer interpolator = InterpolatorType::New();
      resampler->SetInterpolator(interpolator);
      resampler->SetDefaultPixelValue(0);

      IndexImageType::SpacingType spacing;
      spacing[0] = voxelSize;
      spacing[1] = voxelSize;
      spacing[2] = voxelSize;
      resampler->SetOutputSpacing(spacing);

      // Use the same origin
      resampler->SetOutputOrigin(selectIndexImageFilter->GetOutput()->GetOrigin());

      IndexImageType::SizeType    inputSize = baseImage->GetLargestPossibleRegion().GetSize();
      IndexImageType::SpacingType inputSpacing = baseImage->GetSpacing();
      using SizeValueType = IndexImageType::SizeType::SizeValueType;
      IndexImageType::SizeType size;
      size[0] = static_cast<SizeValueType>(inputSize[0] * inputSpacing[0] / voxelSize);
      size[1] = static_cast<SizeValueType>(inputSize[1] * inputSpacing[1] / voxelSize);
      size[2] = static_cast<SizeValueType>(inputSize[2] * inputSpacing[2] / voxelSize);
      resampler->SetSize(size);

      using TransformType = itk::IdentityTransform<double, 3>;
      TransformType::Pointer transform = TransformType::New();
      transform->SetIdentity();
      resampler->SetTransform(transform);
      resampler->Update();
      bvalueImage = resampler->GetOutput();
    }
    else
    {
      bvalueImage = baseImage;
    }
    char        tmpStr[64];
    std::string NrrdValue;
    sprintf(tmpStr, "DWMRI_gradient_%04u", i);
    itk::ExposeMetaData<std::string>(vectorImageReader->GetOutput()->GetMetaDataDictionary(), tmpStr, NrrdValue);
    char tokTmStr[64];
    strcpy(tokTmStr, NrrdValue.c_str());
    TVector tmpDir(3);
    tmpDir[0] = std::stod(strtok(tokTmStr, " "));
    tmpDir[1] = std::stod(strtok(nullptr, " "));
    tmpDir[2] = std::stod(strtok(nullptr, " "));
    if (applyMeasurementFrame)
    {
      std::cout << "Original Direction: " << tmpDir << std::endl;
      tmpDir = measurementFrame * tmpDir;
      std::cout << "New Direction with Measurement Frame: " << tmpDir << std::endl;
    }
    vnl_vector_fixed<double, 3> gradientDir;
    gradientDir[0] = tmpDir[0];
    gradientDir[1] = tmpDir[1];
    gradientDir[2] = tmpDir[2];

    bool useIndex = true;
    for (int j : ignoreIndex)
    {
      if (j == static_cast<int>(i))
      {
        useIndex = false;
      }
    }

    if (useIndex)
    {
      indexImageToVectorImageFilter->SetInput(vectorIndex, bvalueImage);
      gradientDirectionContainer->CreateIndex(vectorIndex);
      gradientDirectionContainer->SetElement(vectorIndex, gradientDir);
      std::cout << "Add Gradient Direction " << vectorIndex << ":  " << gradientDir[0] << ",  " << gradientDir[1]
                << ",  " << gradientDir[2] << std::endl;
      ++vectorIndex;
    }
  }
  indexImageToVectorImageFilter->Update();

  TensorFilterType::Pointer tensorFilter = TensorFilterType::New();
  tensorFilter->SetGradientImage(gradientDirectionContainer, indexImageToVectorImageFilter->GetOutput());
  tensorFilter->SetThreshold(backgroundSuppressingThreshold);
  tensorFilter->SetBValue(BValue);       /* Required */
  tensorFilter->SetNumberOfWorkUnits(1); /* Required */
  if (maskImage.IsNotNull())
  {
    tensorFilter->SetMaskImage(maskImage);
  }
  tensorFilter->Update();

  /* Update the Meta data Header */
  itk::MetaDataDictionary newMeta = tensorFilter->GetOutput()->GetMetaDataDictionary();
  itk::MetaDataDictionary origMeta = vectorImageReader->GetOutput()->GetMetaDataDictionary();
  std::string             NrrdValue;

  itk::ExposeMetaData<std::string>(origMeta, "DWMRI_b-value", NrrdValue);
  itk::EncapsulateMetaData<std::string>(newMeta, "DWMRI_b-value", NrrdValue);

  NrrdValue = "DWMRI";
  itk::EncapsulateMetaData<std::string>(newMeta, "modality", NrrdValue);
  for (int i = 0; i < 4; i++)
  {
    char tmpStr[64];
    sprintf(tmpStr, "NRRD_centerings[%d]", i);
    itk::ExposeMetaData<std::string>(origMeta, tmpStr, NrrdValue);
    itk::EncapsulateMetaData<std::string>(newMeta, tmpStr, NrrdValue);
    sprintf(tmpStr, "NRRD_kinds[%d]", i);
    itk::ExposeMetaData<std::string>(origMeta, tmpStr, NrrdValue);
    itk::EncapsulateMetaData<std::string>(newMeta, tmpStr, NrrdValue);
    sprintf(tmpStr, "NRRD_space units[%d]", i);
    itk::ExposeMetaData<std::string>(origMeta, tmpStr, NrrdValue);
    itk::EncapsulateMetaData<std::string>(newMeta, tmpStr, NrrdValue);
  }

  tensorFilter->GetOutput()->SetMetaDataDictionary(newMeta);

  using TensorImageType = TensorFilterType::TensorImageType;
  TensorImageType::Pointer tensorImage = tensorFilter->GetOutput();

  using WriterType = itk::ImageFileWriter<TensorImageType>;
  WriterType::Pointer nrrdWriter = WriterType::New();
  nrrdWriter->UseCompressionOn();
  nrrdWriter->SetInput(tensorImage);
  nrrdWriter->SetFileName(outputVolume);
  try
  {
    nrrdWriter->Update();
  }
  catch (const itk::ExceptionObject & e)
  {
    std::cout << e << std::endl;
  }
  return EXIT_SUCCESS;
}
