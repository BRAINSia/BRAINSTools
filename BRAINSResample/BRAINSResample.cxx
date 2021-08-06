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
/* ==================================================================
 *
 *  INFO:  NEED TO COMMENT WHAT THIS PROGRAM IS TO BE USED FOR
 *  HACK:  Need to update documentation and licensing.
 *
 *  ================================================================== */

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
#include "BRAINSResampleCLP.h"
#include "BRAINSThreadControl.h"

#include "itkBinaryThresholdImageFilter.h"
#include "itkSignedMaurerDistanceMapImageFilter.h"
#include "itkStatisticsImageFilter.h"
#include "itkMaximumImageFilter.h"
#include "itkDisplacementFieldTransform.h"
#include "GenericTransformImage.h"

#include "TransformToDisplacementField.h"
#include "itkGridForwardWarpImageFilterNew.h"
#include "itkBSplineKernelFunction.h"

#include "itkGridImageSource.h"

#include "BRAINSCommonLib.h"


//// A filter to debug the min/max values
// template <typename TImage>
// void
// PrintImageMinAndMax(TImage * inputImage)
//{
//  using StatisticsFilterType = typename itk::StatisticsImageFilter<TImage>;
//  typename StatisticsFilterType::Pointer statsFilter = StatisticsFilterType::New();
//  statsFilter->SetInput(inputImage);
//  statsFilter->Update();
//  std::cerr << "StatisticsFilter gave Minimum of " << statsFilter->GetMinimum() << " and Maximum of "
//            << statsFilter->GetMaximum() << std::endl;
//}

std::string
inputImageToStringPixelType(const std::string & inputVolume)
{
  // https://itk.org/ITKExamples/src/IO/ImageBase/ReadUnknownImageType/Documentation.html
  itk::ImageIOBase::Pointer imageIO =
    itk::ImageIOFactory::CreateImageIO(inputVolume.c_str(), itk::ImageIOFactory::FileModeType::ReadMode);

  imageIO->SetFileName(inputVolume);
  imageIO->ReadImageInformation();
  using IOPixelType = itk::ImageIOBase::IOPixelType;
  const IOPixelType pixelType = imageIO->GetPixelType();
  std::cout << "Pixel Type is " << itk::ImageIOBase::GetPixelTypeAsString(pixelType) << std::endl;

  using IOComponentType = itk::ImageIOBase::IOComponentEnum;
  const IOComponentType componentType = imageIO->GetComponentType();
  std::cout << "Component Type is " << imageIO->GetComponentTypeAsString(componentType) << std::endl;
  const unsigned int imageDimension = imageIO->GetNumberOfDimensions();
  std::cout << "Image Dimension is " << imageDimension << std::endl;
  if (pixelType == itk::ImageIOBase::SCALAR)
  {
    constexpr int VDimension = 3;
    if (imageDimension == VDimension)
    {
      switch (componentType)
      {
        case itk::ImageIOBase::UCHAR:
        {
          //          using PixelType = unsigned char;
          //          using ImageType = itk::Image<PixelType, VDimension>;
          //          using InternalPixelType = itk::NumericTraits<PixelType>::RealType;
          //          using TBRAINSResampleInternalImageType = itk::Image<InternalPixelType, VDimension>;
          return "uchar";
        }
        case itk::ImageIOBase::CHAR:
        {
          //          using PixelType = char;
          //          using ImageType = itk::Image<PixelType, VDimension>;
          //          using InternalPixelType = itk::NumericTraits<PixelType>::RealType;
          //          using TBRAINSResampleInternalImageType = itk::Image<InternalPixelType, VDimension>;
          return "char";
        }
        case itk::ImageIOBase::USHORT:
        {
          //          using PixelType = unsigned short;
          //          using ImageType = itk::Image<PixelType, VDimension>;
          //          using InternalPixelType = itk::NumericTraits<PixelType>::RealType;
          //          using TBRAINSResampleInternalImageType = itk::Image<InternalPixelType, VDimension>;
          return "ushort";
        }
        case itk::ImageIOBase::SHORT:
        {
          //          using PixelType = short;
          //          using ImageType = itk::Image<PixelType, VDimension>;
          //          using InternalPixelType = itk::NumericTraits<PixelType>::RealType;
          //          using TBRAINSResampleInternalImageType = itk::Image<InternalPixelType, VDimension>;
          return "short";
        }
        case itk::ImageIOBase::UINT:
        {
          //          using PixelType = unsigned int;
          //          using ImageType = itk::Image<PixelType, VDimension>;
          //          using InternalPixelType = itk::NumericTraits<PixelType>::RealType;
          //          using TBRAINSResampleInternalImageType = itk::Image<InternalPixelType, VDimension>;
          return "uint";
        }
        case itk::ImageIOBase::INT:
        {
          //          using PixelType = int;
          //          using ImageType = itk::Image<PixelType, VDimension>;
          //          using InternalPixelType = itk::NumericTraits<PixelType>::RealType;
          //          using TBRAINSResampleInternalImageType = itk::Image<InternalPixelType, VDimension>;
          return "int";
        }
        case itk::ImageIOBase::ULONG:
        {
          //          using PixelType = unsigned long;
          //          using ImageType = itk::Image<PixelType, VDimension>;
          //          using InternalPixelType = itk::NumericTraits<PixelType>::RealType;
          //          using TBRAINSResampleInternalImageType = itk::Image<InternalPixelType, VDimension>;
          return "ulong";
        }
        case itk::ImageIOBase::LONG:
        {
          //          using PixelType = long;
          //          using ImageType = itk::Image<PixelType, VDimension>;
          //          using InternalPixelType = itk::NumericTraits<PixelType>::RealType;
          //          using TBRAINSResampleInternalImageType = itk::Image<InternalPixelType, VDimension>;
          return "long";
        }
        case itk::ImageIOBase::FLOAT:
        {
          //          using PixelType = float;
          //          using ImageType = itk::Image<PixelType, VDimension>;
          //          using InternalPixelType = itk::NumericTraits<PixelType>::RealType;
          //          using TBRAINSResampleInternalImageType = itk::Image<InternalPixelType, VDimension>;
          return "float";
        }
        case itk::ImageIOBase::DOUBLE:
        {
          //          using PixelType = double;
          //          using ImageType = itk::Image<PixelType, VDimension>;
          //          using InternalPixelType = itk::NumericTraits<PixelType>::RealType;
          //          using TBRAINSResampleInternalImageType = itk::Image<InternalPixelType, VDimension>;
          return "double";
        }
        case itk::ImageIOBase::UNKNOWNCOMPONENTTYPE:
          ITK_FALLTHROUGH;
        default:
          std::cerr << "Unknown and unsupported component type!" << std::endl;
          break;
      }
    }
  }
  return "unknown";
}


using InternalPixelType = double; // DEBUG:  This should be "RealType of input pixel type."
using TBRAINSResampleInternalImageType = itk::Image<InternalPixelType, 3>;
using TBRAINSResampleReferenceImageType = TBRAINSResampleInternalImageType;

int
main(int argc, char * argv[])
{
  PARSE_ARGS;
  BRAINSRegisterAlternateIO();
  itk::Object::SetGlobalWarningDisplay(false); // itk warnings aren't thread safe and in
  // this program cause intermittent crashes.

  const BRAINSUtils::StackPushITKDefaultNumberOfThreads TempDefaultNumberOfThreadsHolder(numberOfThreads);

  std::string pixelDataStorageType = "unknown";
  bool        isBinaryImage{ false };
  if (pixelType == "input")
  {
    // Interrogate input image to set the pixelType
    pixelDataStorageType = inputImageToStringPixelType(inputVolume);
  }
  else if (pixelType == "binary")
  {
    pixelDataStorageType = "binary";
    isBinaryImage = true;
  }
  else
  {
    pixelDataStorageType = pixelType;
  }
  pixelType = "PRUPOSEFULLY DESTROYED TYPE TO CAUSE FAILURES"; // pixelType should no longer be used
  // TODO:
  // TODO:  Everything below this line should be wrapped in a "DoIt" template function similar to BRAINSAlignMSP.cxx.

  const bool debug = true;
  bool       useTransform = (!warpTransform.empty());
  const bool useDisplacementField = (!deformationVolume.empty());

  if (inputVolume.empty())
  {
    std::cout << "ERROR: missing input volume name" << std::endl;
    return EXIT_FAILURE;
  }
  if (outputVolume.empty())
  {
    std::cout << "ERROR: missing output volume name" << std::endl;
    return EXIT_FAILURE;
  }
  if (useTransform && useDisplacementField)
  {
    std::cout << "ERROR: warpTransform and deformationVolume are mutually exclusive, only use one of them."
              << std::endl;
    return EXIT_FAILURE;
  }
  // If neither warpTransform nor deformationVolume are defined,
  //  use an identity transform as the warpTransform.
  if (!useTransform && !useDisplacementField)
  {
    std::cout
      << "WARNING: neither warpTransform nor deformationVolume are defined, so warpTransform is set as identity."
      << std::endl;
    useTransform = true;
    warpTransform = "Identity";
  }

  if (debug)
  {
    std::cout << "=====================================================" << std::endl;
    std::cout << "Input Volume:     " << inputVolume << std::endl;
    std::cout << "Reference Volume: " << referenceVolume << std::endl;
    std::cout << "Output Volume:    " << outputVolume << std::endl;
    std::cout << "Pixel Type:       " << pixelDataStorageType << std::endl;
    std::cout << "Interpolation:    " << interpolationMode << std::endl;
    std::cout << "Background Value: " << defaultValue << std::endl;
    if (useDisplacementField)
    {
      std::cout << "Warp by Displacement Volume: " << deformationVolume << std::endl;
    }
    if (useTransform)
    {
      std::cout << "Warp By Transform: " << warpTransform << std::endl;
    }
    std::cout << "=====================================================" << std::endl;
  }
  try
  {
    TBRAINSResampleInternalImageType::Pointer PrincipalOperandImage; // image to be warped
    using ReaderType = itk::ImageFileReader<TBRAINSResampleInternalImageType>;
    ReaderType::Pointer imageReader = ReaderType::New();
    imageReader->SetFileName(inputVolume);
    imageReader->Update();
    PrincipalOperandImage = imageReader->GetOutput();

    // Read ReferenceVolume and DeformationVolume
    using VectorComponentType = double;
    using VectorPixelType = itk::Vector<VectorComponentType, 3>;
    using DisplacementFieldType = itk::Image<VectorPixelType, 3>;
    using DisplacementFieldTransformType = itk::DisplacementFieldTransform<VectorComponentType, 3>;
    // An empty SmartPointer constructor sets up someImage.IsNull() to represent a not-supplied state:
    TBRAINSResampleReferenceImageType::Pointer ReferenceImage;

    ReaderType::Pointer refImageReader = ReaderType::New();
    if (!referenceVolume.empty())
    {
      refImageReader->SetFileName(referenceVolume);
    }
    else
    {
      std::cout << "Warning:  missing Reference Volume defaulted to inputVolume" << std::endl;
      refImageReader->SetFileName(inputVolume);
    }
    refImageReader->Update();
    ReferenceImage = refImageReader->GetOutput();

    // An empty SmartPointer constructor sets up someTransform.IsNull() to
    // represent a not-supplied state:
    itk::Transform<double, 3, 3>::Pointer genericTransform;

    if (useDisplacementField) // it's a warp deformation field
    {
      DisplacementFieldType::Pointer DisplacementField;

      using DefFieldReaderType = itk::ImageFileReader<DisplacementFieldType>;
      DefFieldReaderType::Pointer fieldImageReader = DefFieldReaderType::New();
      fieldImageReader->SetFileName(deformationVolume);
      fieldImageReader->Update();
      DisplacementField = fieldImageReader->GetOutput();
      DisplacementFieldTransformType::Pointer dispTransform = DisplacementFieldTransformType::New();
      dispTransform->SetDisplacementField(DisplacementField.GetPointer());
      genericTransform = dispTransform.GetPointer();
    }
    else if (useTransform)
    {
      try
      {
        if (warpTransform == "Identity")
        {
          itk::VersorRigid3DTransform<double>::Pointer rigidIdentity = itk::VersorRigid3DTransform<double>::New();
          rigidIdentity->SetIdentity();
          genericTransform = rigidIdentity.GetPointer();
        }
        else
        {
          genericTransform = itk::ReadTransformFromDisk(warpTransform);
        }
      }
      catch (itk::ExceptionObject & excp)
      {
        std::cout << "******* HERE *******" << __FILE__ << " " << __LINE__ << excp << std::endl;
        throw excp;
      }
      if (inverseTransform)
      {
        std::string transformFileType = genericTransform->GetNameOfClass();
        std::cout << "Transform File Type:: " << transformFileType << std::endl;
        if (transformFileType == "AffineTransform")
        {
          using LocalAffineTransformType = itk::AffineTransform<double, 3>;
          const LocalAffineTransformType::ConstPointer affineTransform =
            static_cast<LocalAffineTransformType const *>(genericTransform.GetPointer());

          LocalAffineTransformType::Pointer Local_inverseTransform = LocalAffineTransformType::New();
          affineTransform->GetInverse(Local_inverseTransform);

          genericTransform = Local_inverseTransform;
          if (genericTransform.IsNull())
          {
            std::cout << "Error in type conversion " << __FILE__ << __LINE__ << std::endl;
            return EXIT_FAILURE;
          }
        }
        else if (transformFileType == "VersorRigid3DTransform")
        {
          using RigidTransformType = itk::VersorRigid3DTransform<double>;
          const RigidTransformType::ConstPointer rigidTransform =
            static_cast<RigidTransformType const *>(genericTransform.GetPointer());

          RigidTransformType::Pointer Local_inverseTransform = RigidTransformType::New();
          rigidTransform->GetInverse(Local_inverseTransform);

          genericTransform = Local_inverseTransform;

          if (genericTransform.IsNull())
          {
            std::cout << "Error in type conversion " << __FILE__ << __LINE__ << std::endl;
            return EXIT_FAILURE;
          }
        }

        else if (transformFileType == "Similarity3DTransform")
        {
          using Similarity3DTransformType = itk::Similarity3DTransform<double>;
          const Similarity3DTransformType::ConstPointer similarityTransform =
            static_cast<Similarity3DTransformType const *>(genericTransform.GetPointer());

          Similarity3DTransformType::Pointer Local_inverseTransform = Similarity3DTransformType::New();
          similarityTransform->GetInverse(Local_inverseTransform);
          genericTransform = Local_inverseTransform;
          if (genericTransform.IsNull())
          {
            std::cout << "Error in type conversion " << __FILE__ << __LINE__ << std::endl;
            return EXIT_FAILURE;
          }
        }
        else
        {
          std::cout << "*** ERROR ***" << std::endl
                    << " The transform type of " << transformFileType << " does NOT support inverse transformation"
                    << " from BRAINSResample tool." << std::endl;
        }
      }
    }

    TBRAINSResampleInternalImageType::Pointer TransformedImage =
      GenericTransformImage<TBRAINSResampleInternalImageType>(PrincipalOperandImage,
                                                              ReferenceImage,
                                                              // DisplacementField,
                                                              genericTransform.GetPointer(),
                                                              defaultValue,
                                                              interpolationMode,
                                                              isBinaryImage);
    if (gridSpacing.size() == TBRAINSResampleInternalImageType::ImageDimension)
    {
      // find min/max pixels for image
      using StatisticsFilterType = itk::StatisticsImageFilter<TBRAINSResampleInternalImageType>;

      StatisticsFilterType::Pointer statsFilter = StatisticsFilterType::New();
      statsFilter->SetInput(TransformedImage);
      statsFilter->Update();
      TBRAINSResampleInternalImageType::PixelType minPixel(statsFilter->GetMinimum());
      TBRAINSResampleInternalImageType::PixelType maxPixel(statsFilter->GetMaximum());

      DisplacementFieldType::Pointer DisplacementField;
      // create the grid
      if (useTransform)
      { // HACK:  Need to make handeling of transforms more elegant as is done
        // in BRAINSFitHelper.
        using ConverterType = itk::TransformToDisplacementFieldFilter<DisplacementFieldType, double>;
        ConverterType::Pointer myConverter = ConverterType::New();
        myConverter->SetTransform(genericTransform);
        myConverter->SetReferenceImage(TransformedImage);
        myConverter->SetUseReferenceImage(true);
        myConverter->Update();
        DisplacementField = myConverter->GetOutput();
      }
      using MaxFilterType = itk::MaximumImageFilter<TBRAINSResampleInternalImageType>;
      using GFType = itk::GridForwardWarpImageFilterNew<DisplacementFieldType, TBRAINSResampleInternalImageType>;
      GFType::Pointer GFFilter = GFType::New();
      GFFilter->SetInput(DisplacementField);
      GFType::GridSpacingType GridOffsets;
      GridOffsets[0] = gridSpacing[0];
      GridOffsets[1] = gridSpacing[1];
      GridOffsets[2] = gridSpacing[2];
      GFFilter->SetGridPixelSpacing(GridOffsets);
      GFFilter->SetBackgroundValue(minPixel);
      GFFilter->SetForegroundValue(maxPixel);
      // merge grid with warped image
      MaxFilterType::Pointer MFilter = MaxFilterType::New();
      MFilter->SetInput1(GFFilter->GetOutput());
      MFilter->SetInput2(TransformedImage);
      MFilter->Update();
      TransformedImage = MFilter->GetOutput();
    }

    // Write out the output image;  threshold it if necessary.
    if (isBinaryImage)
    {
      // A special case for dealing with binary images
      // where signed distance maps are warped and thresholds created
      using MaskPixelType = short int;
      using MaskImageType = itk::Image<MaskPixelType, 3>;
      using CastImageFilter = itk::CastImageFilter<TBRAINSResampleInternalImageType, MaskImageType>;
      CastImageFilter::Pointer castFilter = CastImageFilter::New();
      castFilter->SetInput(TransformedImage);
      castFilter->Update();

      MaskImageType::Pointer outputImage = castFilter->GetOutput();
      using WriterType = itk::ImageFileWriter<MaskImageType>;
      WriterType::Pointer imageWriter = WriterType::New();
      imageWriter->UseCompressionOn();
      imageWriter->SetFileName(outputVolume);
      imageWriter->SetInput(castFilter->GetOutput());
      try
      {
        imageWriter->Update();
      }
      catch (itk::ExceptionObject & excp)
      {
        std::cout << "******* HERE *******" << __FILE__ << " " << __LINE__ << std::endl;
        std::cout << excp << std::endl;
        return EXIT_FAILURE;
      }
    }
    else if (pixelDataStorageType == "uchar")
    {
      using NewPixelType = unsigned char;
      using NewImageType = itk::Image<NewPixelType, 3>;
      using CastImageFilter = itk::CastImageFilter<TBRAINSResampleInternalImageType, NewImageType>;
      CastImageFilter::Pointer castFilter = CastImageFilter::New();
      castFilter->SetInput(TransformedImage);
      castFilter->Update();

      using WriterType = itk::ImageFileWriter<NewImageType>;
      WriterType::Pointer imageWriter = WriterType::New();
      imageWriter->UseCompressionOn();
      imageWriter->SetFileName(outputVolume);
      imageWriter->SetInput(castFilter->GetOutput());
      try
      {
        imageWriter->Update();
      }
      catch (itk::ExceptionObject & excp)
      {
        std::cout << "******* HERE *******" << __FILE__ << " " << __LINE__ << std::endl;
        std::cout << excp << std::endl;
        return EXIT_FAILURE;
      }
    }
    else if (pixelDataStorageType == "short")
    {
      using NewPixelType = signed short;
      using NewImageType = itk::Image<NewPixelType, 3>;
      using CastImageFilter = itk::CastImageFilter<TBRAINSResampleInternalImageType, NewImageType>;
      CastImageFilter::Pointer castFilter = CastImageFilter::New();
      castFilter->SetInput(TransformedImage);
      castFilter->Update();

      using WriterType = itk::ImageFileWriter<NewImageType>;
      WriterType::Pointer imageWriter = WriterType::New();
      imageWriter->UseCompressionOn();
      imageWriter->SetFileName(outputVolume);
      imageWriter->SetInput(castFilter->GetOutput());
      try
      {
        imageWriter->Update();
      }
      catch (itk::ExceptionObject & excp)
      {
        std::cout << "******* HERE *******" << __FILE__ << " " << __LINE__ << std::endl;
        std::cout << excp << std::endl;
        return EXIT_FAILURE;
      }
    }
    else if (pixelDataStorageType == "ushort")
    {
      using NewPixelType = unsigned short;
      using NewImageType = itk::Image<NewPixelType, 3>;
      using CastImageFilter = itk::CastImageFilter<TBRAINSResampleInternalImageType, NewImageType>;
      CastImageFilter::Pointer castFilter = CastImageFilter::New();
      castFilter->SetInput(TransformedImage);
      castFilter->Update();

      using WriterType = itk::ImageFileWriter<NewImageType>;
      WriterType::Pointer imageWriter = WriterType::New();
      imageWriter->UseCompressionOn();
      imageWriter->SetFileName(outputVolume);
      imageWriter->SetInput(castFilter->GetOutput());
      try
      {
        imageWriter->Update();
      }
      catch (itk::ExceptionObject & excp)
      {
        std::cout << "******* HERE *******" << __FILE__ << " " << __LINE__ << std::endl;
        std::cout << excp << std::endl;
        return EXIT_FAILURE;
      }
    }
    else if (pixelDataStorageType == "int")
    {
      using NewPixelType = int;
      using NewImageType = itk::Image<NewPixelType, 3>;
      using CastImageFilter = itk::CastImageFilter<TBRAINSResampleInternalImageType, NewImageType>;
      CastImageFilter::Pointer castFilter = CastImageFilter::New();
      castFilter->SetInput(TransformedImage);
      castFilter->Update();

      using WriterType = itk::ImageFileWriter<NewImageType>;
      WriterType::Pointer imageWriter = WriterType::New();
      imageWriter->UseCompressionOn();
      imageWriter->SetFileName(outputVolume);
      imageWriter->SetInput(castFilter->GetOutput());
      try
      {
        imageWriter->Update();
      }
      catch (itk::ExceptionObject & excp)
      {
        std::cout << "******* HERE *******" << __FILE__ << " " << __LINE__ << std::endl;
        std::cout << excp << std::endl;
        return EXIT_FAILURE;
      }
    }
    else if (pixelDataStorageType == "uint")
    {
      using NewPixelType = unsigned int;
      using NewImageType = itk::Image<NewPixelType, 3>;
      using CastImageFilter = itk::CastImageFilter<TBRAINSResampleInternalImageType, NewImageType>;
      CastImageFilter::Pointer castFilter = CastImageFilter::New();
      castFilter->SetInput(TransformedImage);
      castFilter->Update();
      using WriterType = itk::ImageFileWriter<NewImageType>;
      WriterType::Pointer imageWriter = WriterType::New();
      imageWriter->UseCompressionOn();
      imageWriter->SetFileName(outputVolume);
      imageWriter->SetInput(castFilter->GetOutput());
      try
      {
        imageWriter->Update();
      }
      catch (itk::ExceptionObject & excp)
      {
        std::cout << "******* HERE *******" << __FILE__ << " " << __LINE__ << std::endl;
        std::cout << excp << std::endl;
        return EXIT_FAILURE;
      }
    }
    else if (pixelDataStorageType == "float")
    {
      using WriterType = itk::ImageFileWriter<TBRAINSResampleInternalImageType>;
      WriterType::Pointer imageWriter = WriterType::New();
      imageWriter->UseCompressionOn();
      imageWriter->SetFileName(outputVolume);
      imageWriter->SetInput(TransformedImage);
      try
      {
        imageWriter->Update();
      }
      catch (itk::ExceptionObject & excp)
      {
        std::cout << "******* HERE *******" << __FILE__ << " " << __LINE__ << std::endl;
        std::cout << excp << std::endl;
        return EXIT_FAILURE;
      }
    }
    else
    {
      std::cout << "ERROR:  Invalid pixelType" << std::endl;
      return EXIT_FAILURE;
    }
  }
  catch (itk::ExceptionObject & excp)
  {
    std::cout << "******* HERE *******" << __FILE__ << " " << __LINE__ << std::endl;
    std::cout << excp << std::endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
