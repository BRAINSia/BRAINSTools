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
#ifndef __itkIO_h
#define __itkIO_h

#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImage.h"
#include "itkCastImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkOrientImageFilter.h"
#include "itkSpatialOrientation.h"
#include "itkSpatialOrientationAdapter.h"
#include "itkMetaDataObject.h"
#include "itkImageRegionIterator.h"
#include "itkThinPlateR2LogRSplineKernelTransform.h"
#include "itkResampleImageFilter.h"
#include "itkImageDuplicator.h"
#include "Imgmath.h"
#include "itkGDCMSeriesFileNames.h"
#include "itkImageSeriesReader.h"
#include "itkGDCMImageIO.h"

namespace itkUtil
{
using SOAdapterType = itk::SpatialOrientationAdapter;
using DirectionType = SOAdapterType::DirectionType;

/**
 *
 *
 *
 */
/** read an image using ITK -- image-based template */
template <typename TImage>
typename TImage::Pointer
ReadImage(const std::string & fileName)
{
  typename TImage::Pointer  image;
  std::string               extension = itksys::SystemTools::GetFilenameLastExtension(fileName);
  itk::GDCMImageIO::Pointer dicomIO = itk::GDCMImageIO::New();
  if (dicomIO->CanReadFile(fileName.c_str()) || (itksys::SystemTools::LowerCase(extension) == ".dcm"))
  {
    std::string dicomDir = itksys::SystemTools::GetParentDirectory(fileName.c_str());

    itk::GDCMSeriesFileNames::Pointer FileNameGenerator = itk::GDCMSeriesFileNames::New();
    FileNameGenerator->SetUseSeriesDetails(true);
    FileNameGenerator->SetDirectory(dicomDir);
    using ContainerType = const std::vector<std::string>;
    const ContainerType & seriesUIDs = FileNameGenerator->GetSeriesUIDs();

    using ReaderType = typename itk::ImageSeriesReader<TImage>;
    typename ReaderType::Pointer reader = ReaderType::New();
    reader->SetFileNames(FileNameGenerator->GetFileNames(seriesUIDs[0]));
    reader->SetImageIO(dicomIO);
    try
    {
      reader->Update();
    }
    catch (itk::ExceptionObject & err)
    {
      std::cout << "Caught an exception: " << std::endl;
      std::cout << err << " " << __FILE__ << " " << __LINE__ << std::endl;
      throw;
    }
    catch (...)
    {
      std::cout << "Error while reading in image for patient " << fileName << std::endl;
      throw;
    }
    image = reader->GetOutput();
    // image->DisconnectPipeline();
    // reader->ReleaseDataFlagOn();
  }
  else
  {
    using ReaderType = itk::ImageFileReader<TImage>;
    typename ReaderType::Pointer reader = ReaderType::New();
    reader->SetFileName(fileName.c_str());
    try
    {
      reader->Update();
    }
    catch (itk::ExceptionObject & err)
    {
      std::cout << "Caught an exception: " << std::endl;
      std::cout << err << " " << __FILE__ << " " << __LINE__ << std::endl;
      throw;
    }
    catch (...)
    {
      std::cout << "Error while reading in image" << fileName << std::endl;
      throw;
    }
    image = reader->GetOutput();
    // image->DisconnectPipeline();
    // reader->ReleaseDataFlagOn();
  }
  return image;
}

/**
 *
 *
 *
 */

template <typename ImageType1, typename ImageType2>
bool
ImagePhysicalDimensionsAreIdentical(typename ImageType1::Pointer & inputImage1,
                                    typename ImageType2::Pointer & inputImage2)
{
  bool same = true;

  same &= (inputImage1->GetDirection() == inputImage2->GetDirection());
  same &= (inputImage1->GetSpacing() == inputImage2->GetSpacing());
  same &= (inputImage1->GetOrigin() == inputImage2->GetOrigin());
  return same;
}

template <typename ImageType>
typename ImageType::Pointer
OrientImage(typename ImageType::ConstPointer &                       inputImage,
            itk::SpatialOrientation::ValidCoordinateOrientationFlags orient)
{
  typename itk::OrientImageFilter<ImageType, ImageType>::Pointer orienter =
    itk::OrientImageFilter<ImageType, ImageType>::New();

  orienter->SetDesiredCoordinateOrientation(orient);
  orienter->UseImageDirectionOn();
  orienter->SetInput(inputImage);
  orienter->Update();
  typename ImageType::Pointer returnval = orienter->GetOutput();
  // returnval->DisconnectPipeline();
  orienter->ReleaseDataFlagOn();
  return returnval;
}

template <typename ImageType>
typename ImageType::Pointer
OrientImage(typename ImageType::ConstPointer & inputImage, const typename ImageType::DirectionType & dirCosines)
{
  return OrientImage<ImageType>(inputImage, SOAdapterType().FromDirectionCosines(dirCosines));
}

template <typename ImageType>
typename ImageType::Pointer
OrientImage(typename ImageType::Pointer & inputImage, const typename ImageType::DirectionType & dirCosines)
{
  typename ImageType::ConstPointer constImg(inputImage);
  return OrientImage<ImageType>(constImg, SOAdapterType().FromDirectionCosines(dirCosines));
}

template <typename ImageType>
typename ImageType::Pointer
OrientImage(typename ImageType::Pointer & inputImage, itk::SpatialOrientation::ValidCoordinateOrientationFlags orient)
{
  typename ImageType::ConstPointer constImg(inputImage);
  return OrientImage<ImageType>(constImg, orient);
}

template <typename ImageType>
typename ImageType::Pointer
ReadImageAndOrient(const std::string & filename, itk::SpatialOrientation::ValidCoordinateOrientationFlags orient)
{
  typename ImageType::Pointer      img = ReadImage<ImageType>(filename);
  typename ImageType::ConstPointer constImg(img);
  typename ImageType::Pointer      image = itkUtil::OrientImage<ImageType>(constImg, orient);
  return image;
}

template <typename ImageType>
typename ImageType::Pointer
ReadImageAndOrient(const std::string & filename, const DirectionType & dir)
{
  return ReadImageAndOrient<ImageType>(filename, SOAdapterType().FromDirectionCosines(dir));
}

template <typename TReadImageType>
typename TReadImageType::Pointer
ReadImageCoronal(const std::string & fileName)
{
  DirectionType CORdir = SOAdapterType().ToDirectionCosines(itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_RIP);

  return ReadImageAndOrient<TReadImageType>(fileName, CORdir);
}

/*
 * This is the prefered ABI for Writing images.
 * We know that the image is not going to change
 * so make sure that the API indicates that.
 */
template <typename ImageType>
void
WriteConstImage(const typename ImageType::ConstPointer image, const std::string & filename)
{
  using WriterType = itk::ImageFileWriter<ImageType>;
  typename WriterType::Pointer writer = WriterType::New();
  writer->UseCompressionOn();
  writer->SetFileName(filename.c_str());
  writer->SetInput(image);
  try
  {
    writer->Update();
  }
  catch (itk::ExceptionObject & err)
  {
    std::cout << "Exception Object caught: " << std::endl;
    std::cout << err << std::endl;
    throw;
  }
}

/*
 * This is for backwards compatibility, and just
 * delegates to WriteConstImage
 * with more restricted write access
 * as necessary.
 */
template <typename ImageType>
void
WriteImage(const typename ImageType::Pointer image, const std::string & filename)
{
  const typename ImageType::ConstPointer temp(image.GetPointer());
  WriteConstImage<ImageType>(temp, filename);
}

/**
 *
 *
 * @author hjohnson (6/4/2008)
 *
 * @param InputImageType
 * @param OutputImageType
 * @param input
 *
 * @return typename OutputImageType::Pointer
 */
template <typename InputImageType, typename OutputImageType>
typename OutputImageType::Pointer
TypeCast(const typename InputImageType::Pointer & input)
{
  using CastToRealFilterType = itk::CastImageFilter<InputImageType, OutputImageType>;
  typename CastToRealFilterType::Pointer toReal = CastToRealFilterType::New();
  toReal->SetInput(input);
  toReal->Update();
  return toReal->GetOutput();
}

/**
 *   \author Hans J. Johnson
 *   Converts images from one type to
 *  another with explicit min and max values. NOTE:  The original
 *  range of the image is determined explicitly from the data,
 *  and then linearly scaled into the range specified.
 * \param image --The input image to convert and scale
 * \param OuputMin --The required minimum value of the output
 *        image
 * \param OutputMax -- The required maximum value of the output
 *        image
 * \return A new image of the specified type and scale.
 */
template <typename InputImageType, typename OutputImageType>
typename OutputImageType::Pointer
ScaleAndCast(const typename InputImageType::Pointer &  image,
             const typename OutputImageType::PixelType OutputMin,
             const typename OutputImageType::PixelType OutputMax)
{
  using R2CRescaleFilterType = itk::RescaleIntensityImageFilter<InputImageType, OutputImageType>;
  typename R2CRescaleFilterType::Pointer RealToProbMapCast = R2CRescaleFilterType::New();
  RealToProbMapCast->SetOutputMinimum(OutputMin);
  RealToProbMapCast->SetOutputMaximum(OutputMax);
  RealToProbMapCast->SetInput(image);
  try
  {
    RealToProbMapCast->Update();
  }
  catch (itk::ExceptionObject & e)
  {
    throw e;
  }
  typename OutputImageType::Pointer returnScaled = RealToProbMapCast->GetOutput();
  return returnScaled;
}

/**
 * This function will do a type cast if the OutputImageType
 * intensity range is larger than the input image type range.
 * If the OutputImageType range is smaller, then a Scaling will
 * occur.
 *
 * @author hjohnson (6/4/2008)
 *
 * @param InputImageType
 * @param OutputImageType
 * @param image
 *
 * @return typename OutputImageType::Pointer
 */
template <typename InputImageType, typename OutputImageType>
typename OutputImageType::Pointer
PreserveCast(const typename InputImageType::Pointer image)
{
  const typename InputImageType::PixelType  inputmin = itk::NumericTraits<typename InputImageType::PixelType>::min();
  const typename InputImageType::PixelType  inputmax = itk::NumericTraits<typename InputImageType::PixelType>::max();
  const typename OutputImageType::PixelType outputmin = itk::NumericTraits<typename OutputImageType::PixelType>::min();
  const typename OutputImageType::PixelType outputmax = itk::NumericTraits<typename OutputImageType::PixelType>::max();
  if ((inputmin >= outputmin) && (inputmax <= outputmax))
  {
    return TypeCast<InputImageType, OutputImageType>(image);
  }
  else
  {
    return ScaleAndCast<InputImageType, OutputImageType>(image, outputmin, outputmax);
  }
}

template <typename ImageType>
typename ImageType::Pointer
CopyImage(const typename ImageType::Pointer & input)
{
  using ImageDupeType = itk::ImageDuplicator<ImageType>;
  typename ImageDupeType::Pointer MyDuplicator = ImageDupeType::New();
  MyDuplicator->SetInputImage(input);
  MyDuplicator->Update();
  return MyDuplicator->GetOutput();
}

/** Common code for allocating an image, allowing the region and spacing to be
 * explicitly set.
 */
template <typename TemplateImageType, typename OutputImageType>
typename OutputImageType::Pointer
AllocateImageFromRegionAndSpacing(const typename TemplateImageType::RegionType &  region,
                                  const typename TemplateImageType::SpacingType & spacing)
{
  typename OutputImageType::Pointer rval;
  rval = OutputImageType::New();
  rval->SetSpacing(spacing);
  //    rval->SetLargestPossibleRegion(region);
  //    rval->SetBufferedRegion(region);
  rval->SetRegions(region);
  rval->Allocate();
  return rval;
}

template <typename ImageType>
typename ImageType::Pointer
AllocateImageFromRegionAndSpacing(const typename ImageType::RegionType &  region,
                                  const typename ImageType::SpacingType & spacing)
{
  return AllocateImageFromRegionAndSpacing<ImageType, ImageType>(region, spacing);
}

/** AllocateImageFromExample creates and allocates an image of the type OutputImageType,
 * using TemplateImageType as the source of size and spacing...
 *
 */
template <typename TemplateImageType, typename OutputImageType>
typename OutputImageType::Pointer
AllocateImageFromExample(const typename TemplateImageType::Pointer & TemplateImage)
{
  typename OutputImageType::Pointer rval = OutputImageType::New();
  rval->CopyInformation(TemplateImage);
  rval->SetRegions(TemplateImage->GetLargestPossibleRegion());
  rval->Allocate();
  return rval;
}

//
// convenience function where template and output images type are the same
template <typename ImageType>
typename ImageType::Pointer
AllocateImageFromExample(const typename ImageType::Pointer & TemplateImage)
{
  return AllocateImageFromExample<ImageType, ImageType>(TemplateImage);
}
} // namespace itkUtil

#endif
