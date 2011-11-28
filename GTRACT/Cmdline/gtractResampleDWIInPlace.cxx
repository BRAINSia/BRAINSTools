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
#include <itkVectorImage.h>
#include <itkImageFileWriter.h>
#include <itkImageFileReader.h>
#include <itkExceptionObject.h>
#include <itkVectorImageToImageAdaptor.h>
#include <itkVectorIndexSelectionCastImageFilter.h>
#include <itkImageRegionIterator.h>
#include <itkImageRegionConstIterator.h>
#include <itkTransformFileWriter.h>

#include "GenericTransformImage.h"
#include "BRAINSFitHelper.h"
#include "BRAINSThreadControl.h"
// #include "itkVectorImageRegisterVersorRigidFilter.h"
// #include "itkVectorImageRegisterAffineFilter.h"

#include "gtractResampleDWIInPlaceCLP.h"
#include "itkImageDuplicator.h"

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
typename IOImageType::Pointer SetVectorImageRigidTransformInPlace(
  typename VersorRigid3DTransformType::ConstPointer RigidTransform, // typename
                                                                    // IOImageType::ConstPointer
                                                                    // InputImage)
  const IOImageType *InputImage)
{
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
  typename IOImageType::Pointer OutputAlignedImage = duplicator->GetOutput();
  // Now change the Origin and Direction to make data aligned.
  OutputAlignedImage->SetOrigin(
    InvOfRigidTransform->GetMatrix() * InputImage->GetOrigin() + InvOfRigidTransform->GetTranslation() );
  OutputAlignedImage->SetDirection( InvOfRigidTransform->GetMatrix() * InputImage->GetDirection() );
  OutputAlignedImage->SetMetaDataDictionary(InputImage->GetMetaDataDictionary() );
  return OutputAlignedImage;
}

int main(int argc, char *argv[])
{
  PARSE_ARGS;
  const BRAINSUtils::StackPushITKDefaultNumberOfThreads TempDefaultNumberOfThreadsHolder(numberOfThreads);

  itk::AddExtraTransformRegister();

  bool debug = true;
  if( debug )
    {
    std::cout << "=====================================================" << std::endl;
    std::cout << "DWI Image: " <<  inputVolume << std::endl;
    std::cout << "Input Transform: " <<  inputTransform << std::endl;
    std::cout << "Output Image: " <<  outputVolume << std::endl;
    std::cout << "Debug Level: " << debugLevel << std::endl;
    std::cout << "=====================================================" << std::endl;
    }

  bool violated = false;
  if( inputVolume.size() == 0 )
    {
    violated = true; std::cout << "  --inputVolume Required! "  << std::endl;
    }
  if( inputTransform.size() == 0 )
    {
    violated = true; std::cout << "  --inputTransform Required! "  << std::endl;
    }
  if( outputVolume.size() == 0 )
    {
    violated = true; std::cout << "  --outputVolume Required! "  << std::endl;
    }
  if( violated )
    {
    return EXIT_FAILURE;
    }

  typedef signed short                        PixelType;
  typedef itk::VectorImage<PixelType, 3>      NrrdImageType;
  typedef itk::VersorRigid3DTransform<double> RigidTransformType;

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
  itk::MetaDataDictionary       inputMetaDataDictionary = resampleImage->GetMetaDataDictionary();
  GenericTransformType::Pointer baseTransform = itk::ReadTransformFromDisk(inputTransform);
  RigidTransformType::Pointer   rigidTransform = dynamic_cast<RigidTransformType *>( baseTransform.GetPointer() );
  if( rigidTransform.IsNull() )
    {
    std::cerr << "Transform "
              << inputTransform
              << " wasn't of the expected type itk::VersorRigid3DTransform<double"
              << std::endl;
    return EXIT_FAILURE;
    }
  resampleImage = SetVectorImageRigidTransformInPlace<NrrdImageType>(rigidTransform.GetPointer(), resampleImage);
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

    //    std::cout << &curGradientDirection[0] << " " <<
    // &curGradientDirection[1] << " " << &curGradientDirection[2] << std::endl;

    // Rotated the Diffusion Gradient
    curGradientDirection = rigidTransform->GetRotationMatrix().GetVnlMatrix() * curGradientDirection;

    // Updated the Image MetaData Dictionary with Updated Gradient Information
    sprintf(tmpStr, " %18.15lf %18.15lf %18.15lf", curGradientDirection[0], curGradientDirection[1],
            curGradientDirection[2]);
    NrrdValue = tmpStr;
    itk::EncapsulateMetaData<std::string>(resampleImage->GetMetaDataDictionary(), KeyString, NrrdValue);
    }

  typedef itk::ImageFileWriter<NrrdImageType> WriterType;
  WriterType::Pointer nrrdWriter = WriterType::New();
  nrrdWriter->UseCompressionOn();
  nrrdWriter->UseInputMetaDataDictionaryOn();
  nrrdWriter->SetInput( resampleImage );
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
