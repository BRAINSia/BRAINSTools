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
/*==================================================================

INFO:  NEED TO COMMENT WHAT THIS PROGRAM IS TO BE USED FOR

==================================================================*/

#include <iostream>
#include "itkVector.h"
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkWarpImageFilter.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkNearestNeighborInterpolateImageFunction.h"
#include "itkBSplineInterpolateImageFunction.h"
#include "itkWindowedSincInterpolateImageFunction.h"
#include "ApplyWarpCLP.h"

#include "itkBinaryThresholdImageFilter.h"
#include "itkSignedMaurerDistanceMapImageFilter.h"
#include "itkStatisticsImageFilter.h"
#include "itkIO.h"

#include "GenericTransformImage.h"

// A filter to debug the min/max values
template <typename TImage>
void
PrintImageMinAndMax(TImage * inputImage)
{
  //  typename TImage::PixelType resultMaximum:
  //  typename TImage::PixelType resultMinimum;
  using StatisticsFilterType = typename itk::StatisticsImageFilter<TImage>;
  typename StatisticsFilterType::Pointer statsFilter = StatisticsFilterType::New();
  statsFilter->SetInput(inputImage);
  statsFilter->Update();
  //  resultMaximum = statsFilter->GetMaximum();
  //  resultMinimum = statsFilter->GetMinimum();
  std::cerr << "StatisticsFilter gave Minimum of " << statsFilter->GetMinimum() << " and Maximum of "
            << statsFilter->GetMaximum() << std::endl;
}

int
ApplyWarp(int argc, char * argv[])
{
  PARSE_ARGS;
  BRAINSRegisterAlternateIO();

  const bool useTransform = (warpTransform.size() > 0);
  {
    const bool useDisplacementField = (deformationVolume.size() > 0);
    const bool debug = true;

    if (debug)
    {
      std::cout << "=====================================================" << std::endl;
      std::cout << "Input Volume:      " << inputVolume << std::endl;
      std::cout << "Reference Volume:  " << referenceVolume << std::endl;
      std::cout << "Output Volume:     " << outputVolume << std::endl;
      std::cout << "Pixel Type:        " << pixelType << std::endl;
      std::cout << "Orientation to RAI:" << orientationRAI << std::endl;
      std::cout << "Interpolation:     " << interpolationMode << std::endl;
      std::cout << "Background Value:  " << defaultValue << std::endl;
      if (useDisplacementField)
      {
        std::cout << "Warp by Deformation Volume: " << deformationVolume << std::endl;
      }
      if (useTransform)
      {
        std::cout << "Warp By Transform: " << warpTransform << std::endl;
      }
      std::cout << "=====================================================" << std::endl;
    }

    if (useTransformMode.size() > 0)
    {
      std::cout << "Scripting 'code rot' note:  The useTransformMode parameter will be ignored.  Now ApplyWarp infers "
                   "the warpTransform type from the contents of the .mat file."
                << std::endl;
    }

    if (useTransform == useDisplacementField)
    {
      std::cout
        << "Choose one of the two possibilities, a BRAINSFit transform --or-- a high-dimensional deformation field."
        << std::endl;
      exit(1);
    }
  }

  using ImageType = itk::Image<float, 3>;
  using RefImageType = itk::Image<float, 3>;
  ImageType::Pointer PrincipalOperandImage; // One name for the image to be warped.
  {

    if (orientationRAI)
    {
      PrincipalOperandImage = itkUtil::ReadImage<ImageType>(inputVolume);
      PrincipalOperandImage =
        itkUtil::OrientImage<ImageType>(PrincipalOperandImage, itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_RAI);
    }
    else
    {
      PrincipalOperandImage = itkUtil::ReadImage<ImageType>(inputVolume);
    }
  }

  // Read ReferenceVolume and DeformationVolume

  using VectorComponentType = float;
  using VectorPixelType = itk::Vector<VectorComponentType, 3>;
  using DisplacementFieldType = itk::Image<VectorPixelType, 3>;

  // An empty SmartPointer constructor sets up someImage.IsNull() to represent a not-supplied state:
  DisplacementFieldType::Pointer DisplacementField;
  RefImageType::Pointer          ReferenceImage;

  if (useTransform)
  {

    if (referenceVolume.size() > 0)
    {
      if (orientationRAI)
      {
        ReferenceImage = itkUtil::ReadImage<RefImageType>(referenceVolume);
        ReferenceImage =
          itkUtil::OrientImage<RefImageType>(ReferenceImage, itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_RAI);
      }
      else
      {
        ReferenceImage = itkUtil::ReadImage<RefImageType>(referenceVolume);
      }
    }
    else
    {
      std::cout << "Alert:  missing Reference Volume defaulted to: " << inputVolume << std::endl;
      //  ReferenceImage = itkUtil::ReadImage<RefImageType>( inputVolume );
      if (orientationRAI)
      {
        ReferenceImage = itkUtil::ReadImage<RefImageType>(inputVolume);
        ReferenceImage =
          itkUtil::OrientImage<RefImageType>(ReferenceImage, itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_RAI);
      }
      else
      {
        ReferenceImage = itkUtil::ReadImage<RefImageType>(inputVolume);
      }
    }
  }
  else if (!useTransform) // that is, it's a warp by deformation field:
  {
    if (orientationRAI)
    {
      DisplacementField = itkUtil::ReadImage<DisplacementFieldType>(deformationVolume);
      DisplacementField = itkUtil::OrientImage<DisplacementFieldType>(
        DisplacementField, itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_RAI);
    }
    else
    {
      DisplacementField = itkUtil::ReadImage<DisplacementFieldType>(deformationVolume);
    }
    if (referenceVolume.size() > 0)
    {
      if (orientationRAI)
      {
        ReferenceImage = itkUtil::ReadImage<RefImageType>(referenceVolume);
        ReferenceImage =
          itkUtil::OrientImage<RefImageType>(ReferenceImage, itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_RAI);
      }
      else
      {
        ReferenceImage = itkUtil::ReadImage<RefImageType>(referenceVolume);
      }
    }
  }

  // Read optional transform:

  // An empty SmartPointer constructor sets up someTransform.IsNull() to represent a not-supplied state:
  BSplineTransformType::Pointer itkBSplineTransform;
  AffineTransformType::Pointer  ITKAffineTransform;

  if (useTransform)
  {
    std::cerr << "Invalid option, not implemented in ITKv4 version yet." << std::endl;
    return -1;
  }

  ImageType::Pointer TransformedImage =
    GenericTransformImage<ImageType, RefImageType, DisplacementFieldType>(PrincipalOperandImage,
                                                                          ReferenceImage,
                                                                          DisplacementField,
                                                                          NULL,
                                                                          defaultValue,
                                                                          interpolationMode,
                                                                          pixelType == "binary");

  // Write out the output image;  threshold it if necessary.
  if (pixelType == "binary")
  {
    // A special case for dealing with binary images
    // where signed distance maps are warped and thresholds created
    using MaskPixelType = short int;
    using MaskImageType = itk::Image<MaskPixelType, 3>;
    using CastImageFilter = itk::CastImageFilter<ImageType, MaskImageType>;
    CastImageFilter::Pointer castFilter = CastImageFilter::New();
    castFilter->SetInput(TransformedImage);
    castFilter->Update();

    MaskImageType::Pointer outputImage = castFilter->GetOutput();
    using WriterType = itk::ImageFileWriter<MaskImageType>;
    WriterType::Pointer imageWriter = WriterType::New();
    imageWriter->SetFileName(outputVolume);
    imageWriter->SetInput(castFilter->GetOutput());
    imageWriter->Update();
  }
  else if (pixelType == "uchar")
  {
    using NewPixelType = unsigned char;
    using NewImageType = itk::Image<NewPixelType, 3>;
    using CastImageFilter = itk::CastImageFilter<ImageType, NewImageType>;
    CastImageFilter::Pointer castFilter = CastImageFilter::New();
    castFilter->SetInput(TransformedImage);
    castFilter->Update();

    using WriterType = itk::ImageFileWriter<NewImageType>;
    WriterType::Pointer imageWriter = WriterType::New();
    imageWriter->SetFileName(outputVolume);
    imageWriter->SetInput(castFilter->GetOutput());
    imageWriter->Update();
  }
  else if (pixelType == "short")
  {
    using NewPixelType = signed short;
    using NewImageType = itk::Image<NewPixelType, 3>;
    using CastImageFilter = itk::CastImageFilter<ImageType, NewImageType>;
    CastImageFilter::Pointer castFilter = CastImageFilter::New();
    castFilter->SetInput(TransformedImage);
    castFilter->Update();

    using WriterType = itk::ImageFileWriter<NewImageType>;
    WriterType::Pointer imageWriter = WriterType::New();
    imageWriter->SetFileName(outputVolume);
    imageWriter->SetInput(castFilter->GetOutput());
    imageWriter->Update();
  }
  else if (pixelType == "ushort")
  {
    using NewPixelType = unsigned short;
    using NewImageType = itk::Image<NewPixelType, 3>;
    using CastImageFilter = itk::CastImageFilter<ImageType, NewImageType>;
    CastImageFilter::Pointer castFilter = CastImageFilter::New();
    castFilter->SetInput(TransformedImage);
    castFilter->Update();

    using WriterType = itk::ImageFileWriter<NewImageType>;
    WriterType::Pointer imageWriter = WriterType::New();
    imageWriter->SetFileName(outputVolume);
    imageWriter->SetInput(castFilter->GetOutput());
    imageWriter->Update();
  }
  else if (pixelType == "int")
  {
    using NewPixelType = int;
    using NewImageType = itk::Image<NewPixelType, 3>;
    using CastImageFilter = itk::CastImageFilter<ImageType, NewImageType>;
    CastImageFilter::Pointer castFilter = CastImageFilter::New();
    castFilter->SetInput(TransformedImage);
    castFilter->Update();

    using WriterType = itk::ImageFileWriter<NewImageType>;
    WriterType::Pointer imageWriter = WriterType::New();
    imageWriter->SetFileName(outputVolume);
    imageWriter->SetInput(castFilter->GetOutput());
    imageWriter->Update();
  }
  else if (pixelType == "uint")
  {
    using NewPixelType = unsigned int;
    using NewImageType = itk::Image<NewPixelType, 3>;
    using CastImageFilter = itk::CastImageFilter<ImageType, NewImageType>;
    CastImageFilter::Pointer castFilter = CastImageFilter::New();
    castFilter->SetInput(TransformedImage);
    castFilter->Update();
    ;
    using WriterType = itk::ImageFileWriter<NewImageType>;
    WriterType::Pointer imageWriter = WriterType::New();
    imageWriter->SetFileName(outputVolume);
    imageWriter->SetInput(castFilter->GetOutput());
    imageWriter->Update();
  }
  else if (pixelType == "float")
  {
    using WriterType = itk::ImageFileWriter<ImageType>;
    WriterType::Pointer imageWriter = WriterType::New();
    imageWriter->SetFileName(outputVolume);
    imageWriter->SetInput(TransformedImage);
    imageWriter->Update();
  }
  else
  {
    std::cout << "ERROR:  Invalid pixelType" << std::endl;
    exit(-1);
  }

  return EXIT_SUCCESS;
}

int
main(int argc, char * argv[])
{
  std::cout << "This program has been replaced by BRAINSResample.  PLEASE TRY TO AVOID USING THIS!" << std::endl;
  std::cout << "This program has been replaced by BRAINSResample.  PLEASE TRY TO AVOID USING THIS!" << std::endl;
  std::cout << "This program has been replaced by BRAINSResample.  PLEASE TRY TO AVOID USING THIS!" << std::endl;
  std::cout << "This program has been replaced by BRAINSResample.  PLEASE TRY TO AVOID USING THIS!" << std::endl;
  std::cout << "This program has been replaced by BRAINSResample.  PLEASE TRY TO AVOID USING THIS!" << std::endl;
  std::cout << "This program has been replaced by BRAINSResample.  PLEASE TRY TO AVOID USING THIS!" << std::endl;
  std::cout << "This program has been replaced by BRAINSResample.  PLEASE TRY TO AVOID USING THIS!" << std::endl;
  std::cout << "This program has been replaced by BRAINSResample.  PLEASE TRY TO AVOID USING THIS!" << std::endl;

  return ApplyWarp(argc, argv);
}
