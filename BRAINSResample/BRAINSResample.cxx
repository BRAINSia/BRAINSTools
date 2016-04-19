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
 *  TODO:  NEED TO COMMENT WHAT THIS PROGRAM IS TO BE USED FOR
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

typedef float                                                                  InternalPixelType;
typedef itk::Image<InternalPixelType, 3> TBRAINSResampleInternalImageType;
typedef TBRAINSResampleInternalImageType                                       TBRAINSResampleReferenceImageType;

// A filter to debug the min/max values
template <class TImage>
void PrintImageMinAndMax(TImage *inputImage)
{
  typedef typename itk::StatisticsImageFilter<TImage> StatisticsFilterType;
  typename StatisticsFilterType::Pointer statsFilter = StatisticsFilterType::New();
  statsFilter->SetInput(inputImage);
  statsFilter->Update();
  std::cerr << "StatisticsFilter gave Minimum of " << statsFilter->GetMinimum()
            << " and Maximum of " << statsFilter->GetMaximum() << std::endl;
}

int main(int argc, char *argv[])
{
  PARSE_ARGS;
  BRAINSRegisterAlternateIO();
  itk::Object::SetGlobalWarningDisplay(false); // itk warnings aren't thread safe and in
  // this program cause intermittent crashes.

  const BRAINSUtils::StackPushITKDefaultNumberOfThreads TempDefaultNumberOfThreadsHolder(numberOfThreads);

  const bool debug = true;
  bool useTransform = ( warpTransform.size() > 0 );
  const bool useDisplacementField = ( deformationVolume.size() > 0 );

  if(inputVolume.empty())
    {
    std::cout << "ERROR: missing input volume name"
              << std::endl;
    return EXIT_FAILURE;
    }
  if(outputVolume.empty())
    {
    std::cout << "ERROR: missing output volume name"
              << std::endl;
    return EXIT_FAILURE;
    }
  if( useTransform && useDisplacementField )
    {
    std::cout << "ERROR: warpTransform and deformationVolume are mutually exclusive, only use one of them."
              << std::endl;
    return EXIT_FAILURE;
    }
  // If neither warpTransform nor deformationVolume are defined,
  //  use an identity transform as the warpTransform.
  if( !useTransform && !useDisplacementField )
    {
    std::cout << "WARNING: neither warpTransform nor deformationVolume are defined, so warpTransform is set as identity."
              << std::endl;
    useTransform = true;
    warpTransform = "Identity";
    }

  if( debug )
    {
    std::cout << "=====================================================" << std::endl;
    std::cout << "Input Volume:     " <<  inputVolume << std::endl;
    std::cout << "Reference Volume: " <<  referenceVolume << std::endl;
    std::cout << "Output Volume:    " <<  outputVolume << std::endl;
    std::cout << "Pixel Type:       " <<  pixelType << std::endl;
    std::cout << "Interpolation:    " <<  interpolationMode << std::endl;
    std::cout << "Background Value: " <<  defaultValue << std::endl;
    if( useDisplacementField )
      {
      std::cout << "Warp by Displacement Volume: " << deformationVolume << std::endl;
      }
    if( useTransform )
      {
      std::cout << "Warp By Transform: " << warpTransform << std::endl;
      }
    std::cout << "=====================================================" << std::endl;
    }
  try
    {
    TBRAINSResampleInternalImageType::Pointer PrincipalOperandImage;  // image to be warped
    typedef itk::ImageFileReader<TBRAINSResampleInternalImageType> ReaderType;
    ReaderType::Pointer imageReader = ReaderType::New();
    imageReader->SetFileName(inputVolume);
    imageReader->Update();
    PrincipalOperandImage = imageReader->GetOutput();

    // Read ReferenceVolume and DeformationVolume
    typedef double                                                                    VectorComponentType;
    typedef itk::Vector<VectorComponentType, 3> VectorPixelType;
    typedef itk::Image<VectorPixelType,  3>     DisplacementFieldType;
    typedef itk::DisplacementFieldTransform<VectorComponentType,3>
      DisplacementFieldTransformType;
    // An empty SmartPointer constructor sets up someImage.IsNull() to represent a not-supplied state:
    TBRAINSResampleReferenceImageType::Pointer ReferenceImage;

    ReaderType::Pointer refImageReader = ReaderType::New();
    if( referenceVolume.size() > 0 )
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

    if( useDisplacementField )  // it's a warp deformation field
      {
      DisplacementFieldType::Pointer DisplacementField;

      typedef itk::ImageFileReader<DisplacementFieldType> DefFieldReaderType;
      DefFieldReaderType::Pointer fieldImageReader = DefFieldReaderType::New();
      fieldImageReader->SetFileName( deformationVolume );
      fieldImageReader->Update();
      DisplacementField = fieldImageReader->GetOutput();
      DisplacementFieldTransformType::Pointer dispTransform =
        DisplacementFieldTransformType::New();
      dispTransform->SetDisplacementField(DisplacementField.GetPointer());
      genericTransform = dispTransform.GetPointer();
      }
    else if( useTransform )
      {
      try
        {
        if ( warpTransform == "Identity" )
          {
          itk::VersorRigid3DTransform<double>::Pointer rigidIdentity = itk::VersorRigid3DTransform<double>::New();
          rigidIdentity->SetIdentity();
          genericTransform = rigidIdentity.GetPointer();
          }
        else
          {
          genericTransform = itk::ReadTransformFromDisk( warpTransform );
          }
        }
      catch( itk::ExceptionObject & excp )
        {
        std::cout << "******* HERE *******" << __FILE__ << " " << __LINE__ << excp << std::endl;
        throw excp;
        }
      if( inverseTransform )
        {
        std::string transformFileType = genericTransform->GetNameOfClass();
        std::cout << "Transform File Type:: " << transformFileType << std::endl;
        if( transformFileType == "AffineTransform" )
          {
          typedef itk::AffineTransform<double, 3>
            LocalAffineTransformType;
          const LocalAffineTransformType::ConstPointer affineTransform =
            static_cast<LocalAffineTransformType const *>(
              genericTransform.GetPointer() );

          LocalAffineTransformType::Pointer Local_inverseTransform = LocalAffineTransformType::New();
          affineTransform->GetInverse( Local_inverseTransform );

          genericTransform = Local_inverseTransform;
          if( genericTransform.IsNull() )
            {
            std::cout << "Error in type conversion " << __FILE__ << __LINE__ << std::endl;
            return EXIT_FAILURE;
            }
          }
        else if( transformFileType == "VersorRigid3DTransform" )
          {
          typedef itk::VersorRigid3DTransform<double>
            RigidTransformType;
          const RigidTransformType::ConstPointer rigidTransform =
            static_cast<RigidTransformType const *>(
              genericTransform.GetPointer() );

          RigidTransformType::Pointer Local_inverseTransform = RigidTransformType::New();
          rigidTransform->GetInverse( Local_inverseTransform );

          genericTransform = Local_inverseTransform;

          if( genericTransform.IsNull() )
            {
            std::cout << "Error in type conversion " << __FILE__ << __LINE__ << std::endl;
            return EXIT_FAILURE;
            }
          }
        else
          {
          std::cout << "*** ERROR ***" << std::endl
            << " The transform type of " << transformFileType
            << " does NOT support inverse transformation" << std::endl;
          }
        }
      }

    TBRAINSResampleInternalImageType::Pointer TransformedImage =
      GenericTransformImage<TBRAINSResampleInternalImageType, TBRAINSResampleInternalImageType, DisplacementFieldType>(
        PrincipalOperandImage,
        ReferenceImage,
        // DisplacementField,
        genericTransform.GetPointer(),
        defaultValue,
        interpolationMode,
        pixelType == "binary");
    if( gridSpacing.size() == TBRAINSResampleInternalImageType::ImageDimension )
      {
      // find min/max pixels for image
      typedef itk::StatisticsImageFilter<TBRAINSResampleInternalImageType> StatisticsFilterType;

      StatisticsFilterType::Pointer statsFilter =
        StatisticsFilterType::New();
      statsFilter->SetInput(TransformedImage);
      statsFilter->Update();
      TBRAINSResampleInternalImageType::PixelType minPixel( statsFilter->GetMinimum() );
      TBRAINSResampleInternalImageType::PixelType maxPixel( statsFilter->GetMaximum() );

      DisplacementFieldType::Pointer DisplacementField;
      // create the grid
      if( useTransform )
        { // HACK:  Need to make handeling of transforms more elegant as is done
          // in BRAINSFitHelper.
        typedef itk::TransformToDisplacementFieldFilter<DisplacementFieldType, double> ConverterType;
        ConverterType::Pointer myConverter = ConverterType::New();
        myConverter->SetTransform(genericTransform);
        myConverter->SetReferenceImage(TransformedImage);
        myConverter->SetUseReferenceImage(true);
        myConverter->Update();
        DisplacementField = myConverter->GetOutput();
        }
      typedef itk::MaximumImageFilter<TBRAINSResampleInternalImageType> MaxFilterType;
      typedef itk::GridForwardWarpImageFilterNew
        <DisplacementFieldType, TBRAINSResampleInternalImageType> GFType;
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
      MFilter->SetInput1( GFFilter->GetOutput() );
      MFilter->SetInput2(TransformedImage);
      MFilter->Update();
      TransformedImage = MFilter->GetOutput();
      }

    // Write out the output image;  threshold it if necessary.
    if( pixelType == "binary" )
      {
      // A special case for dealing with binary images
      // where signed distance maps are warped and thresholds created
      typedef short int                                                             MaskPixelType;
      typedef itk::Image<MaskPixelType,  3>   MaskImageType;
      typedef itk::CastImageFilter<TBRAINSResampleInternalImageType, MaskImageType> CastImageFilter;
      CastImageFilter::Pointer castFilter = CastImageFilter::New();
      castFilter->SetInput(TransformedImage);
      castFilter->Update();

      MaskImageType::Pointer outputImage = castFilter->GetOutput();
      typedef itk::ImageFileWriter<MaskImageType> WriterType;
      WriterType::Pointer imageWriter = WriterType::New();
      imageWriter->UseCompressionOn();
      imageWriter->SetFileName(outputVolume);
      imageWriter->SetInput( castFilter->GetOutput() );
      try
        {
        imageWriter->Update();
        }
      catch( itk::ExceptionObject & excp )
        {
        std::cout << "******* HERE *******" << __FILE__ << " " << __LINE__ << std::endl;
        std::cout << excp << std::endl;
        return EXIT_FAILURE;
        }
      }
    else if( pixelType == "uchar" )
      {
      typedef unsigned char                                                        NewPixelType;
      typedef itk::Image<NewPixelType, 3>    NewImageType;
      typedef itk::CastImageFilter<TBRAINSResampleInternalImageType, NewImageType> CastImageFilter;
      CastImageFilter::Pointer castFilter = CastImageFilter::New();
      castFilter->SetInput(TransformedImage);
      castFilter->Update();

      typedef itk::ImageFileWriter<NewImageType> WriterType;
      WriterType::Pointer imageWriter = WriterType::New();
      imageWriter->UseCompressionOn();
      imageWriter->SetFileName(outputVolume);
      imageWriter->SetInput( castFilter->GetOutput() );
      try
        {
        imageWriter->Update();
        }
      catch( itk::ExceptionObject & excp )
        {
        std::cout << "******* HERE *******" << __FILE__ << " " << __LINE__ << std::endl;
        std::cout << excp << std::endl;
        return EXIT_FAILURE;
        }
      }
    else if( pixelType == "short" )
      {
      typedef signed short                                                         NewPixelType;
      typedef itk::Image<NewPixelType, 3>    NewImageType;
      typedef itk::CastImageFilter<TBRAINSResampleInternalImageType, NewImageType> CastImageFilter;
      CastImageFilter::Pointer castFilter = CastImageFilter::New();
      castFilter->SetInput(TransformedImage);
      castFilter->Update();

      typedef itk::ImageFileWriter<NewImageType> WriterType;
      WriterType::Pointer imageWriter = WriterType::New();
      imageWriter->UseCompressionOn();
      imageWriter->SetFileName(outputVolume);
      imageWriter->SetInput( castFilter->GetOutput() );
      try
        {
        imageWriter->Update();
        }
      catch( itk::ExceptionObject & excp )
        {
        std::cout << "******* HERE *******" << __FILE__ << " " << __LINE__ << std::endl;
        std::cout << excp << std::endl;
        return EXIT_FAILURE;
        }
      }
    else if( pixelType == "ushort" )
      {
      typedef unsigned short                                                       NewPixelType;
      typedef itk::Image<NewPixelType, 3>    NewImageType;
      typedef itk::CastImageFilter<TBRAINSResampleInternalImageType, NewImageType> CastImageFilter;
      CastImageFilter::Pointer castFilter = CastImageFilter::New();
      castFilter->SetInput(TransformedImage);
      castFilter->Update();

      typedef itk::ImageFileWriter<NewImageType> WriterType;
      WriterType::Pointer imageWriter = WriterType::New();
      imageWriter->UseCompressionOn();
      imageWriter->SetFileName(outputVolume);
      imageWriter->SetInput( castFilter->GetOutput() );
      try
        {
        imageWriter->Update();
        }
      catch( itk::ExceptionObject & excp )
        {
        std::cout << "******* HERE *******" << __FILE__ << " " << __LINE__ << std::endl;
        std::cout << excp << std::endl;
        return EXIT_FAILURE;
        }
      }
    else if( pixelType == "int" )
      {
      typedef int                                                                  NewPixelType;
      typedef itk::Image<NewPixelType, 3>    NewImageType;
      typedef itk::CastImageFilter<TBRAINSResampleInternalImageType, NewImageType> CastImageFilter;
      CastImageFilter::Pointer castFilter = CastImageFilter::New();
      castFilter->SetInput(TransformedImage);
      castFilter->Update();

      typedef itk::ImageFileWriter<NewImageType> WriterType;
      WriterType::Pointer imageWriter = WriterType::New();
      imageWriter->UseCompressionOn();
      imageWriter->SetFileName(outputVolume);
      imageWriter->SetInput( castFilter->GetOutput() );
      try
        {
        imageWriter->Update();
        }
      catch( itk::ExceptionObject & excp )
        {
        std::cout << "******* HERE *******" << __FILE__ << " " << __LINE__ << std::endl;
        std::cout << excp << std::endl;
        return EXIT_FAILURE;
        }
      }
    else if( pixelType == "uint" )
      {
      typedef unsigned int                                                         NewPixelType;
      typedef itk::Image<NewPixelType, 3>    NewImageType;
      typedef itk::CastImageFilter<TBRAINSResampleInternalImageType, NewImageType> CastImageFilter;
      CastImageFilter::Pointer castFilter = CastImageFilter::New();
      castFilter->SetInput(TransformedImage);
      castFilter->Update();
      typedef itk::ImageFileWriter<NewImageType> WriterType;
      WriterType::Pointer imageWriter = WriterType::New();
      imageWriter->UseCompressionOn();
      imageWriter->SetFileName(outputVolume);
      imageWriter->SetInput( castFilter->GetOutput() );
      try
        {
        imageWriter->Update();
        }
      catch( itk::ExceptionObject & excp )
        {
        std::cout << "******* HERE *******" << __FILE__ << " " << __LINE__ << std::endl;
        std::cout << excp << std::endl;
        return EXIT_FAILURE;
        }
      }
    else if( pixelType == "float" )
      {
      typedef itk::ImageFileWriter<TBRAINSResampleInternalImageType> WriterType;
      WriterType::Pointer imageWriter = WriterType::New();
      imageWriter->UseCompressionOn();
      imageWriter->SetFileName(outputVolume);
      imageWriter->SetInput(TransformedImage);
      try
        {
        imageWriter->Update();
        }
      catch( itk::ExceptionObject & excp )
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
  catch( itk::ExceptionObject & excp )
    {
    std::cout << "******* HERE *******" << __FILE__ << " " << __LINE__ << std::endl;
    std::cout << excp << std::endl;
    return EXIT_FAILURE;
    }
  return EXIT_SUCCESS;
}
