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
#include <cmath>

#include <itkImage.h>
#include <itkVectorImage.h>
#include <itkImageFileWriter.h>
#include <itkImageFileReader.h>
#include <itkExceptionObject.h>
#include <itkVectorIndexSelectionCastImageFilter.h>
#include <itkComposeImageFilter.h>
#include <itkImageRegionIterator.h>
#include <itkImageRegionConstIterator.h>
#include <itkResampleImageFilter.h>
#include <itkVectorResampleImageFilter.h>
#include <itkTransformFileWriter.h>
#include <itkMetaDataDictionary.h>
#include <itkMetaDataObject.h>

#include "GenericTransformImage.h"
#include "BRAINSFitHelper.h"
#include "BRAINSThreadControl.h"

#include "gtractResampleDWIInPlaceCLP.h"
#include "itkImageDuplicator.h"
#include "itkNumberToString.h"
/**
 * \author Hans J. Johnson
 * \brief This templated function will duplicate the input image, change the direction and origin to refelect the physical space
 * tranform that would be equivalent to calling the resample image filter.
 * InplaceImage=SetVectorImageRigidTransformInPlace(RigidTransform,InputImage); ResampleImage(InplaceImage,Identity);
 * should produce the same result as ResampleImage(InputImage,RigidTransform);
 * \param RigidTransform -- Currently must be a VersorRigid3D
 * \param InputImage The image to be duplicated and modified to incorporate the rigid transform.
 * \return an image with the same voxels values as the input, but with differnt physical space representation.
 */
template <class IOImageType>
typename IOImageType::Pointer
SetVectorImageRigidTransformInPlace(typename itk::VersorRigid3DTransform<double>::ConstPointer RigidTransform,
                                    const IOImageType *InputImage)
{
  typedef itk::VersorRigid3DTransform<double>              VersorRigid3DTransformType;

  typename VersorRigid3DTransformType::Pointer InvOfRigidTransform = VersorRigid3DTransformType::New();
  const typename IOImageType::PointType centerPoint = RigidTransform->GetCenter();
  InvOfRigidTransform->SetCenter( centerPoint );
  InvOfRigidTransform->SetIdentity();
  RigidTransform->GetInverse(InvOfRigidTransform);

  /** Wei: The output image will have exact the same index contents
    but with modified image info so that the index-to-physical mapping
    makes the image in the physical space AC-PC aligned */
  typedef itk::ImageDuplicator<IOImageType> DuplicatorType;
  typename DuplicatorType::Pointer duplicator = DuplicatorType::New();
  duplicator->SetInputImage(InputImage);
  duplicator->Update();
  typename IOImageType::Pointer OutputAlignedImage = duplicator->GetModifiableOutput();
  // Now change the Origin and Direction to make data aligned.
  OutputAlignedImage->SetOrigin(
    InvOfRigidTransform->GetMatrix() * InputImage->GetOrigin() + InvOfRigidTransform->GetOffset() );
  OutputAlignedImage->SetDirection( InvOfRigidTransform->GetMatrix() * InputImage->GetDirection() );
  OutputAlignedImage->SetMetaDataDictionary(InputImage->GetMetaDataDictionary() );
  return OutputAlignedImage;
}

int main(int argc, char *argv[])
{
  PARSE_ARGS;
  BRAINSRegisterAlternateIO();
  const BRAINSUtils::StackPushITKDefaultNumberOfThreads TempDefaultNumberOfThreadsHolder(numberOfThreads);
  itk::NumberToString<double>                                        doubleConvert;
  typedef itk::Transform<double, 3, 3> GenericTransformType;

  bool debug = true;
  if( debug )
    {
    std::cout << "=====================================================" << std::endl;
    std::cout << "DWI Image: " <<  inputVolume << std::endl;
    std::cout << "Input Transform: " <<  inputTransform << std::endl;
    std::cout << "warpDWI Transform: " <<  warpDWITransform << std::endl;
    std::cout << "Output Image: " <<  outputVolume << std::endl;
    std::cout << "Debug Level: " << debugLevel << std::endl;
    std::cout << "Image Output Size: "
              << imageOutputSize[0] << "," << imageOutputSize[1] << "," << imageOutputSize[2] << std::endl;
    std::cout << "=====================================================" << std::endl;
    }

  bool violated = false;
  if( inputVolume.size() == 0 )
    {
    violated = true;
    std::cout << "  --inputVolume Required! "  << std::endl;
    }
  if( inputTransform.size() == 0  && warpDWITransform.size() == 0 )
    {
    violated = true;
    std::cout << "  --inputTransform or --warpDWITransform Required! "  << std::endl;
    }
  if( outputVolume.size() == 0 )
    {
    violated = true;
    std::cout << "  --outputVolume Required! "  << std::endl;
    }
  if( violated )
    {
    return EXIT_FAILURE;
    }

  typedef signed short                        PixelType;
  typedef itk::VectorImage<PixelType, 3>      NrrdImageType;
  typedef itk::VersorRigid3DTransform<double> RigidTransformType;
  typedef itk::Image<PixelType, 3>            SingleComponentImageType;
  typedef itk::ImageFileReader<NrrdImageType,
                               itk::DefaultConvertPixelTraits<PixelType> > FileReaderType;
  FileReaderType::Pointer imageReader = FileReaderType::New();
  imageReader->SetFileName( inputVolume );
  try
    {
    imageReader->Update();
    }
  catch( itk::ExceptionObject & ex )
    {
    std::cout << ex << std::endl;
    throw;
    }
  std::cout << "Read Image" << std::endl;

  NrrdImageType::Pointer        resampleImage = imageReader->GetOutput();
  NrrdImageType::DirectionType  myDirection = resampleImage->GetDirection();
  itk::MetaDataDictionary       inputMetaDataDictionary = resampleImage->GetMetaDataDictionary();
  GenericTransformType::Pointer baseTransform = ITK_NULLPTR;
  if( inputTransform == "ID"  || inputTransform == "Identity" || inputTransform.size() == 0 )
    {
    RigidTransformType::Pointer LocalRigidTransform = RigidTransformType::New();
    LocalRigidTransform->SetIdentity();
    baseTransform = LocalRigidTransform;
    }
  else
    {
    baseTransform = itk::ReadTransformFromDisk(inputTransform);
    }
  RigidTransformType::Pointer rigidTransform = dynamic_cast<RigidTransformType *>( baseTransform.GetPointer() );

  if( rigidTransform.IsNull() )
    {
    std::cerr << "Transform "
              << inputTransform
              << " wasn't of the expected type itk::VersorRigid3DTransform<double>"
              << std::endl;
    return EXIT_FAILURE;
    }

  GenericTransformType::Pointer     warpDWIXFRM = ITK_NULLPTR;
  if(warpDWITransform.size() > 0)
    {
    if ( warpDWITransform.find(".nii") != std::string::npos )
      {
      // Read ReferenceVolume and DeformationVolume
      typedef double                                                VectorComponentType;
      typedef itk::Vector<VectorComponentType, 3> VectorPixelType;
      typedef itk::Image<VectorPixelType,  3>     DisplacementFieldType;
      typedef itk::DisplacementFieldTransform<VectorComponentType,3> DisplacementFieldTransformType;

      typedef itk::ImageFileReader<DisplacementFieldType> DefFieldReaderType;
      DefFieldReaderType::Pointer fieldImageReader = DefFieldReaderType::New();
      fieldImageReader->SetFileName( warpDWITransform );
      fieldImageReader->Update();

      DisplacementFieldType::Pointer DisplacementField = fieldImageReader->GetOutput();

      DisplacementFieldTransformType::Pointer dispTransform =
        DisplacementFieldTransformType::New();
      dispTransform->SetDisplacementField(DisplacementField.GetPointer());
      warpDWIXFRM = dispTransform.GetPointer();
      }
    else
      {
      warpDWIXFRM = itk::ReadTransformFromDisk(warpDWITransform);
      }
    }
  // Get measurement frame and its inverse from DWI scan
  std::vector<std::vector<double> > msrFrame;
  itk::ExposeMetaData<std::vector<std::vector<double> > >( inputMetaDataDictionary,
                                                           "NRRD_measurement frame", msrFrame );
  vnl_matrix_fixed<double, 3, 3> DWIMeasurementFrame;
  if( msrFrame.size() != 0 )
    {
    for( unsigned int i = 0; i < 3; i++ )
      {
      for( unsigned int j = 0; j < 3; j++ )
        {
        DWIMeasurementFrame[i][j] = msrFrame[i][j];
        }
      }
    }
  else
    {
    std::cout << "File does not have NRRD measurement frame metadata!"
              << std::endl
              << "If this is not a DWI NRRD file (e.g. fMRI NIFTI, etc.), you can safely ignore this"
              << std::endl;
    DWIMeasurementFrame.set_identity();
    }
  vnl_matrix_fixed<double, 3, 3> DWIInverseMeasurementFrame = vnl_inverse( DWIMeasurementFrame );

  // Resample DWI in place
  resampleImage = SetVectorImageRigidTransformInPlace<NrrdImageType>(rigidTransform.GetPointer(), resampleImage);

  std::cout << "Rigid transform matrix: " << rigidTransform->GetMatrix().GetVnlMatrix() << std::endl;

  std::stringstream outputGradDirMetaDataStream;

  // Rotate gradient vectors by rigid transform and inverse measurement frame
  for( unsigned int i = 0; i < resampleImage->GetVectorLength(); i++ )
    {
    // Get Current Gradient Direction
    vnl_vector<double> curGradientDirection(3);
    char               tmpStr[64];
    sprintf(tmpStr, "DWMRI_gradient_%04u", i);
    std::string KeyString(tmpStr);
    std::string NrrdValue;

    itk::ExposeMetaData<std::string>(inputMetaDataDictionary, KeyString, NrrdValue);
    sscanf(
      NrrdValue.c_str(), " %lf %lf %lf", &curGradientDirection[0], &curGradientDirection[1], &curGradientDirection[2]);

    // Rotate the diffusion gradient with rigid transform and inverse measurement frame
    RigidTransformType::Pointer inverseRigidTransform = RigidTransformType::New();
    const NrrdImageType::PointType centerPoint = rigidTransform->GetCenter();
    inverseRigidTransform->SetCenter( centerPoint );
    inverseRigidTransform->SetIdentity();
    rigidTransform->GetInverse(inverseRigidTransform);

    curGradientDirection = inverseRigidTransform->GetMatrix().GetVnlMatrix() * DWIInverseMeasurementFrame
      * curGradientDirection;

    // Updated the Image MetaData Dictionary with Updated Gradient Information
    // sprintf(tmpStr, " %18.15lf %18.15lf %18.15lf", curGradientDirection[0], curGradientDirection[1],
    // curGradientDirection[2]);
    // NrrdValue = tmpStr;
    NrrdValue = " ";
    for( unsigned dir = 0; dir < 3; ++dir )
      {
      if( dir > 0 )
        {
        NrrdValue += " ";
        }
      NrrdValue += doubleConvert(curGradientDirection[dir]);
      }
    itk::EncapsulateMetaData<std::string>(resampleImage->GetMetaDataDictionary(), KeyString, NrrdValue);

    outputGradDirMetaDataStream << KeyString << ","
                                << curGradientDirection[0] << ","
                                << curGradientDirection[1] << ","
                                << curGradientDirection[2] << std::endl;
    }

  // Set DWI measurement frame to identity by multiplying by its inverse
  // Update Image MetaData Dictionary with new measurement frame
  vnl_matrix_fixed<double, 3, 3>    newMeasurementFrame = DWIInverseMeasurementFrame * DWIMeasurementFrame;
  std::vector<std::vector<double> > newMf(3);
  std::stringstream outputMFMetaDataStream;
  outputMFMetaDataStream << "measurement frame";
  for( unsigned int i = 0; i < 3; i++ )
    {
    newMf[i].resize(3);
    for( unsigned int j = 0; j < 3; j++ )
      {
      newMf[i][j] = newMeasurementFrame[i][j];
      outputMFMetaDataStream << "," << newMf[i][j];
      }
    }
  outputMFMetaDataStream << std::endl;
  itk::EncapsulateMetaData<std::vector<std::vector<double> > >(
    resampleImage->GetMetaDataDictionary(), "NRRD_measurement frame", newMf );

  // If --writeOutputMetaData ./metaData.csv is specified on the command line,
  // then write out the output image measurement frame and diffusion gradient directions in a simple CSV file.
  // This helps to verify the output image meta data in TestSuite.
  if( writeOutputMetaData != "" )
    {
    std::ofstream outputMetaDataFile;
    outputMetaDataFile.open( writeOutputMetaData.c_str() );
    if( !outputMetaDataFile.is_open() )
      {
      std::cerr << "Can't write the output meta data CSV file " << writeOutputMetaData << std::endl;
      }
    outputMetaDataFile << outputMFMetaDataStream.str();
    outputMetaDataFile << outputGradDirMetaDataStream.str();
    outputMetaDataFile.close();
    }

  // Pad image
  const NrrdImageType::RegionType    inputRegion = resampleImage->GetLargestPossibleRegion();
  const NrrdImageType::SizeType      inputSize = inputRegion.GetSize();
  const NrrdImageType::PointType     inputOrigin = resampleImage->GetOrigin();
  const NrrdImageType::DirectionType inputDirection = resampleImage->GetDirection();
  const NrrdImageType::SpacingType   inputSpacing = resampleImage->GetSpacing();

  NrrdImageType::SizeType newSize;
  std::vector<size_t>     imagePadding(6, 0);
  for( int qq = 0; qq < 3; ++qq )
    {
    if( ( imageOutputSize[qq] > 0 ) && ( static_cast<itk::SizeValueType>( imageOutputSize[qq] )  > inputSize[qq] ) )
      {
      const size_t sizeDiff = imageOutputSize[qq] - inputSize[qq];
      const size_t halfPadding = sizeDiff / 2;
      const size_t isOddDiff = sizeDiff % 2;
      imagePadding[qq] =  halfPadding  + isOddDiff;
      imagePadding[qq + 3] = halfPadding;
      }
    newSize[qq] = inputSize[qq] + imagePadding[qq] + imagePadding[qq + 3];
    }

  vnl_matrix_fixed<double, 3, 3> inputDirectionMatrix;
  vnl_matrix_fixed<double, 3, 3> inputSpacingMatrix;
  for( unsigned int i = 0; i < 3; i++ )
    {
    for( unsigned int j = 0; j < 3; j++ )
      {
      inputDirectionMatrix[i][j] = inputDirection[i][j];
      if( i == j )
        {
        inputSpacingMatrix[i][j] = inputSpacing[i];
        }
      else
        {
        inputSpacingMatrix[i][j] = 0.0;
        }
      }
    }

  vnl_matrix_fixed<double, 3, 3> spaceDirections = inputDirectionMatrix * inputSpacingMatrix;
  vnl_matrix_fixed<double, 3, 4> newMatrix;
  for( unsigned int i = 0; i < 3; i++ )
    {
    for( unsigned int j = 0; j < 4; j++ )
      {
      if( j == 3 )
        {
        newMatrix[i][j] = inputOrigin[i];
        }
      else
        {
        newMatrix[i][j] = spaceDirections[i][j];
        }
      }
    }

  vnl_matrix_fixed<double, 4, 1> voxelShift;
  voxelShift[0][0] = -1.0 * ( double )( imagePadding[0] );
  voxelShift[1][0] = -1.0 * ( double )( imagePadding[1] );
  voxelShift[2][0] = -1.0 * ( double )( imagePadding[2] );
  voxelShift[3][0] = 1.0;

  vnl_matrix_fixed<double, 3, 1> newOriginMatrix = newMatrix * voxelShift;

  NrrdImageType::PointType newOrigin;
  for( unsigned int i = 0; i < 3; i++ )
    {
    newOrigin[i] = newOriginMatrix[i][0];
    }

  std::cout << "Input DWI Image Origin: ( " << inputOrigin[0] << ", " << inputOrigin[1] << ", " << inputOrigin[2]
            << " )" << std::endl;
  std::cout << "Input DWI Image Size: " << inputSize[0] << " " << inputSize[1] << " " << inputSize[2] << std::endl;
  std::cout << " " << std::endl;
  std::cout << "Output DWI Image Origin: ( " << newOrigin[0] << ", " << newOrigin[1] << ", " << newOrigin[2] << " )"
            << std::endl;
  std::cout << "Output DWI Image Size: " << newSize[0] << " " << newSize[1] << " " << newSize[2] << std::endl;

  NrrdImageType::Pointer paddedImage = NrrdImageType::New();
  paddedImage->CopyInformation(resampleImage);
  paddedImage->SetVectorLength( resampleImage->GetVectorLength() );
  paddedImage->SetMetaDataDictionary( resampleImage->GetMetaDataDictionary() );
  paddedImage->SetRegions( newSize );
  paddedImage->SetOrigin( newOrigin );
  paddedImage->Allocate();

  NrrdImageType::Pointer finalImage;
  typedef itk::ImageRegionIterator<NrrdImageType> IteratorType;
  IteratorType InIt( resampleImage, resampleImage->GetRequestedRegion() );
  for( InIt.GoToBegin(); !InIt.IsAtEnd(); ++InIt )
    {
    NrrdImageType::IndexType InIndex = InIt.GetIndex();
    NrrdImageType::IndexType OutIndex;
    OutIndex[0] = InIndex[0] + imagePadding[0];
    OutIndex[1] = InIndex[1] + imagePadding[1];
    OutIndex[2] = InIndex[2] + imagePadding[2];

    NrrdImageType::PixelType InImagePixel = resampleImage->GetPixel( InIndex );
    paddedImage->SetPixel( OutIndex, InImagePixel );
    }
  paddedImage->Update();

  if(referenceVolume != "")
    {
    //For each component, extract, resample to list, and finally compose back into a vector image.
    typedef itk::ImageFileReader<SingleComponentImageType> ReferenceFileReaderType;
    ReferenceFileReaderType::Pointer referenceImageReader = ReferenceFileReaderType::New();
    referenceImageReader->SetFileName(referenceVolume);
    referenceImageReader->Update();

    const size_t lengthOfPixelVector = paddedImage->GetVectorLength();

    typedef itk::ComposeImageFilter<SingleComponentImageType, NrrdImageType> ComposeCovariantVectorImageFilterType;
    ComposeCovariantVectorImageFilterType::Pointer composer= ComposeCovariantVectorImageFilterType::New();

    typedef itk::VectorIndexSelectionCastImageFilter< NrrdImageType, SingleComponentImageType >
        VectorIndexSelectionCastImageFilterType;
    VectorIndexSelectionCastImageFilterType::Pointer vectorImageToImageSelector = VectorIndexSelectionCastImageFilterType::New();
    vectorImageToImageSelector->SetInput(paddedImage);
    for(size_t componentToExtract = 0 ; componentToExtract < lengthOfPixelVector; ++componentToExtract )
      {
      vectorImageToImageSelector->SetIndex( componentToExtract );
      vectorImageToImageSelector->Update();

      // Resample to a new space with basic linear/identity transform.
      typedef itk::ResampleImageFilter<SingleComponentImageType,SingleComponentImageType> ComponentResamplerType;
      ComponentResamplerType::Pointer componentResampler = ComponentResamplerType::New();
      componentResampler->SetOutputParametersFromImage(referenceImageReader->GetOutput());
      componentResampler->SetInput(vectorImageToImageSelector->GetOutput());
      if(warpDWIXFRM.IsNotNull())
        {
        componentResampler->SetTransform(warpDWIXFRM);
        }
      //default to linear
      //default to IdentityTransform
      // TODO feed it composite transform.
      //default background value of 0.
      componentResampler->Update();
      //Add to list for Compose
      composer->SetInput(componentToExtract,componentResampler->GetOutput());
      }
    composer->Update();
    finalImage = composer->GetOutput();
    finalImage->SetMetaDataDictionary(paddedImage->GetMetaDataDictionary() );
    }
  else
    {
    finalImage=paddedImage;
    }

  // Write out resampled in place DWI
  typedef itk::ImageFileWriter<NrrdImageType> WriterType;
  WriterType::Pointer nrrdWriter = WriterType::New();
  nrrdWriter->UseCompressionOn();
  nrrdWriter->UseInputMetaDataDictionaryOn();
  nrrdWriter->SetInput( finalImage );
  nrrdWriter->SetFileName( outputVolume );
  try
    {
    nrrdWriter->Update();
    }
  catch( itk::ExceptionObject & e )
    {
    std::cout << e << std::endl;
    }
  if(!outputResampledB0.empty())
    {
    const size_t B0Index=0;

    typedef itk::VectorIndexSelectionCastImageFilter< NrrdImageType, SingleComponentImageType >
      VectorIndexSelectionCastImageFilterType;
    VectorIndexSelectionCastImageFilterType::Pointer vectorImageToImageSelector = VectorIndexSelectionCastImageFilterType::New();
    vectorImageToImageSelector->SetInput(finalImage);
    vectorImageToImageSelector->SetIndex( B0Index );
    vectorImageToImageSelector->Update();

    // Write out resampled in place DWI
    typedef itk::ImageFileWriter<SingleComponentImageType> B0WriterType;
    B0WriterType::Pointer B0Writer = B0WriterType::New();
    B0Writer->UseCompressionOn();
    B0Writer->SetInput( vectorImageToImageSelector->GetOutput() );
    B0Writer->SetFileName( outputResampledB0 );
    try
      {
      B0Writer->Update();
      }
    catch( itk::ExceptionObject & e )
      {
      std::cout << e << std::endl;
      }
    }
  return EXIT_SUCCESS;
}
