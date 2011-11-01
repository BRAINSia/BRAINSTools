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
#include <itkImageToVectorImageFilter.h>
#include <itkExtractImageFilter.h>
#include <itkImageMaskSpatialObject.h>

#include <itkDiffusionTensor3DReconstructionWithMaskImageFilter.h>
#include "itkGtractImageIO.h"
#include "itkGtractParameterIO.h"
#include "itkComputeDiffusionTensorImageFilter.h"
#include "itkVectorResampleImageFilter.h"
#include "itkIdentityTransform.h"
#include "itkMetaDataObject.h"
#include "itkIOCommon.h"
#include "gtractTensorCLP.h"
#include "BRAINSThreadControl.h"

#include "GenericTransformImage.h"
#include "itkBRAINSROIAutoImageFilter.h"
#include "itkIO.h"

int main(int argc, char *argv[])
{
  PARSE_ARGS;
  BRAINSUtils::SetThreadCount(numberOfThreads);
  itk::AddExtraTransformRegister();

  typedef signed short                   PixelType;
  typedef double                         TensorPixelType;
  typedef itk::VectorImage<PixelType, 3> VectorImageType;
  typedef itk::Image<PixelType, 3>       IndexImageType;
  typedef itk::Image<unsigned char, 3>   MaskImageType;
  IndexImageType::SizeType MedianFilterSize;
  MedianFilterSize[0] = medianFilterSize[0];
  MedianFilterSize[1] = medianFilterSize[1];
  MedianFilterSize[2] = medianFilterSize[2];

  bool debug = true;
  // applyMeasurementFrame = true;
  if( debug )
    {
    std::cout << "=====================================================" << std::endl;
    std::cout << "Input Image: " <<  inputVolume << std::endl;
    std::cout << "Output Image: " <<  outputVolume << std::endl;
    std::cout << "Resample Isotropic: " << resampleIsotropic << std::endl;
    std::cout << "Voxel Size: " << voxelSize << std::endl;
    std::cout << "Median Filter Size: " << MedianFilterSize << std::endl;
    std::cout << "Threshold: " << backgroundSuppressingThreshold << std::endl;
    std::cout << "B0 Index: " << b0Index << std::endl;
    std::cout << "Apply Measurement Frame: " << applyMeasurementFrame << std::endl;
    std::cout << "=====================================================" << std::endl;
    }

  bool violated = false;
  if( inputVolume.size() == 0 )
    {
    violated = true; std::cout << "  --inputVolume Required! "  << std::endl;
    }
  if( outputVolume.size() == 0 )
    {
    violated = true; std::cout << "  --outputVolume Required! "  << std::endl;
    }
  if( violated )
    {
    return EXIT_FAILURE;
    }

  typedef itk::ImageFileReader<VectorImageType,
                               itk::DefaultConvertPixelTraits<PixelType> > VectorImageReaderType;
  VectorImageReaderType::Pointer vectorImageReader = VectorImageReaderType::New();
  vectorImageReader->SetFileName( inputVolume );

  try
    {
    vectorImageReader->Update();
    }
  catch( itk::ExceptionObject & ex )
    {
    std::cout << ex << std::endl;
    throw;
    }
  // figure out mask processing
  MaskImageType::ConstPointer maskImage; // will stay NULL if no mask is used.

  if( maskProcessingMode == "ROIAUTO" )
    {
    typedef itk::Image<PixelType, 3> IndexImageType;
    typedef itk::VectorIndexSelectionCastImageFilter
      <VectorImageType, IndexImageType> VectorSelectFilterType;
    VectorSelectFilterType::Pointer SelectIndexImageFilter =
      VectorSelectFilterType::New();
    SelectIndexImageFilter->SetIndex(b0Index);
    SelectIndexImageFilter->SetInput(vectorImageReader->GetOutput() );
    try
      {
      SelectIndexImageFilter->Update();
      }
    catch( itk::ExceptionObject e )
      {
      std::cerr << e << std::endl;
      return EXIT_FAILURE;
      }
    ;
    IndexImageType::Pointer b0Image(SelectIndexImageFilter->GetOutput() );
    typedef itk::BRAINSROIAutoImageFilter<IndexImageType, MaskImageType>
      ROIAutoType;
    ROIAutoType::Pointer ROIFilter = ROIAutoType::New();
    ROIFilter->SetInput(b0Image);
    ROIFilter->SetClosingSize(9.0);  // default TODO Make parameter
    ROIFilter->SetDilateSize(0.0);
    ROIFilter->Update();
    maskImage = ROIFilter->GetBinaryImageROI();
    }
  else if( maskProcessingMode == "ROI" )
    {
    if( maskVolume == "" )
      {
      std::cerr << "Error: missing mask Volume needed for ROI mask Processing"
                << std::endl;
      return EXIT_FAILURE;
      }
    maskImage = itkUtil::ReadImage<MaskImageType>(maskVolume);
    if( maskImage.IsNull() )
      {
      std::cerr << "Error: can't read mask volume "
                << maskVolume << std::endl;
      return EXIT_FAILURE;
      }
    }
  /* Extract Diffusion Information from the Header */
  std::string BValue_str;
  std::string BValue_keyStr("DWMRI_b-value");
  itk::ExposeMetaData<std::string>(vectorImageReader->GetOutput()->GetMetaDataDictionary(),
                                   BValue_keyStr.c_str(), BValue_str);
  double BValue = atof( BValue_str.c_str() );
  std::cout << "The BValue was found to be " << BValue_str << std::endl;

  std::vector<std::vector<double> > msrFrame;
  itk::ExposeMetaData<std::vector<std::vector<double> > >(
    vectorImageReader->GetOutput()->GetMetaDataDictionary(),
    "NRRD_measurement frame", msrFrame);
  TMatrix measurementFrame(3, 3);
  for( int i = 0; i < 3; i++ )
    {
    for( int j = 0; j < 3; j++ )
      {
      measurementFrame[i][j] = msrFrame[i][j];
      }
    }

  /* Process Invidual B-value Images and Reassemble the Vector Image */
  typedef itk::DiffusionTensor3DReconstructionWithMaskImageFilter<PixelType, PixelType,
                                                                  TensorPixelType> TensorFilterType;
  typedef TensorFilterType::GradientDirectionContainerType
    DirectionContainerType;
  DirectionContainerType::Pointer gradientDirectionContainer = DirectionContainerType::New();

  typedef itk::ImageToVectorImageFilter<IndexImageType> VectorImageFilterType;
  VectorImageFilterType::Pointer indexImageToVectorImageFilter = VectorImageFilterType::New();
  int                            vectorIndex = 0;
  for( unsigned int i = 0; i < vectorImageReader->GetOutput()->GetVectorLength(); i++ )
    {
    typedef itk::VectorIndexSelectionCastImageFilter<VectorImageType, IndexImageType> VectorSelectFilterType;
    typedef VectorSelectFilterType::Pointer                                           VectorSelectFilterPointer;

    VectorSelectFilterPointer selectIndexImageFilter = VectorSelectFilterType::New();
    selectIndexImageFilter->SetIndex( i );
    selectIndexImageFilter->SetInput( vectorImageReader->GetOutput() );
    try
      {
      selectIndexImageFilter->Update();
      }
    catch( itk::ExceptionObject e )
      {
      std::cout << e << std::endl;
      }

    /* Median Filter */
    IndexImageType::Pointer baseImage;
    if( MedianFilterSize[0] > 0  ||  MedianFilterSize[1] > 0  ||  MedianFilterSize[2] > 0 )
      {
      typedef itk::MedianImageFilter<IndexImageType, IndexImageType> MedianFilterType;
      MedianFilterType::Pointer filter = MedianFilterType::New();
      filter->SetInput( selectIndexImageFilter->GetOutput() );
      filter->SetRadius( MedianFilterSize );
      filter->Update();
      baseImage = filter->GetOutput();
      }
    else
      {
      baseImage = selectIndexImageFilter->GetOutput();
      }

    /* Resample To Isotropic Images */
    IndexImageType::Pointer bvalueImage;
    if( resampleIsotropic )
      {
      typedef itk::ResampleImageFilter<IndexImageType, IndexImageType> ResampleFilterType;
      ResampleFilterType::Pointer resampler = ResampleFilterType::New();
      resampler->SetInput( baseImage );

      typedef itk::LinearInterpolateImageFunction<IndexImageType, double> InterpolatorType;
      InterpolatorType::Pointer interpolator = InterpolatorType::New();
      resampler->SetInterpolator( interpolator );
      resampler->SetDefaultPixelValue( 0 );

      IndexImageType::SpacingType spacing;
      spacing[0] = voxelSize;
      spacing[1] = voxelSize;
      spacing[2] = voxelSize;
      resampler->SetOutputSpacing( spacing );

      // Use the same origin
      resampler->SetOutputOrigin( selectIndexImageFilter->GetOutput()->GetOrigin() );

      IndexImageType::SizeType    inputSize  = baseImage->GetLargestPossibleRegion().GetSize();
      IndexImageType::SpacingType inputSpacing  = baseImage->GetSpacing();
      typedef IndexImageType::SizeType::SizeValueType SizeValueType;
      IndexImageType::SizeType size;
      size[0] = static_cast<SizeValueType>( inputSize[0] * inputSpacing[0] / voxelSize );
      size[1] = static_cast<SizeValueType>( inputSize[1] * inputSpacing[1] / voxelSize );
      size[2] = static_cast<SizeValueType>( inputSize[2] * inputSpacing[2] / voxelSize );
      resampler->SetSize( size );

      typedef itk::IdentityTransform<double, 3> TransformType;
      TransformType::Pointer transform = TransformType::New();
      transform->SetIdentity();
      resampler->SetTransform( transform );
      resampler->Update();
      bvalueImage = resampler->GetOutput();
      }
    else
      {
      bvalueImage = baseImage;
      }
    char        tmpStr[64];
    std::string NrrdValue;
    sprintf(tmpStr, "DWMRI_gradient_%04d", i);
    itk::ExposeMetaData<std::string>(vectorImageReader->GetOutput()->GetMetaDataDictionary(), tmpStr, NrrdValue);
    char tokTmStr[64];
    strcpy( tokTmStr, NrrdValue.c_str() );
    TVector tmpDir(3);
    tmpDir[0] = atof( strtok(tokTmStr, " ") );
    tmpDir[1] = atof( strtok(NULL, " ") );
    tmpDir[2] = atof( strtok(NULL, " ") );
    if( applyMeasurementFrame )
      {
      std::cout << "Original Direction: " << tmpDir << std::endl;
      tmpDir = measurementFrame * tmpDir;
      std::cout << "New Direction with Measurement Frame: " << tmpDir << std::endl;
      }
    vnl_vector_fixed<double, 3> gradientDir;
    gradientDir[0] = tmpDir[0]; gradientDir[1] = tmpDir[1]; gradientDir[2] = tmpDir[2];

    bool useIndex = true;
    for( unsigned int j = 0; j < ignoreIndex.size(); j++ )
      {
      if( ignoreIndex.at(j) == static_cast<int>( i ) )
        {
        useIndex = false;
        }
      }

    if( useIndex )
      {
      indexImageToVectorImageFilter->SetInput( vectorIndex, bvalueImage );
      gradientDirectionContainer->CreateIndex(vectorIndex);
      gradientDirectionContainer->SetElement(vectorIndex, gradientDir);
      std::cout << "Add Gradient Direction " << vectorIndex << ":  " << gradientDir[0] << ",  " << gradientDir[1]
                << ",  " << gradientDir[2] << std::endl;
      vectorIndex++;
      }
    }
  indexImageToVectorImageFilter->Update();

  TensorFilterType::Pointer tensorFilter = TensorFilterType::New();
  tensorFilter->SetGradientImage( gradientDirectionContainer, indexImageToVectorImageFilter->GetOutput() );
  tensorFilter->SetThreshold( backgroundSuppressingThreshold );
  tensorFilter->SetBValue(BValue);     /* Required */
  tensorFilter->SetNumberOfThreads(1); /* Required */
  if( maskImage.IsNotNull() )
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
  for( int i = 0; i < 4; i++ )
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

  typedef TensorFilterType::TensorImageType TensorImageType;
  TensorImageType::Pointer tensorImage = tensorFilter->GetOutput();

  typedef itk::ImageFileWriter<TensorImageType> WriterType;
  WriterType::Pointer nrrdWriter = WriterType::New();
  nrrdWriter->UseCompressionOn();
  nrrdWriter->SetInput( tensorImage );
  nrrdWriter->SetFileName( outputVolume );
  try
    {
    nrrdWriter->Update();
    }
  catch( itk::ExceptionObject e )
    {
    std::cout << e << std::endl;
    }
  return EXIT_SUCCESS;
}
