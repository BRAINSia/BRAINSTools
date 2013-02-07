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
#include <itkVectorImageToImageAdaptor.h>
#include <itkVectorIndexSelectionCastImageFilter.h>
#include <itkImageRegionIterator.h>
#include <itkImageRegionConstIterator.h>
#include <itkTransformFileWriter.h>
#include <itkMetaDataDictionary.h>
#include <itkMetaDataObject.h>
#include "GenericTransformImage.h"
#include "BRAINSFitHelper.h"
#include "BRAINSThreadControl.h"
// #include "itkVectorImageRegisterVersorRigidFilter.h"
// #include "itkVectorImageRegisterAffineFilter.h"

#include "gtractResampleDWIInPlaceCLP.h"
#include "itkImageDuplicator.h"
#include "DoubleToString.h"
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
  DoubleToString                                        doubleConvert;
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
  NrrdImageType::DirectionType  myDirection = resampleImage->GetDirection();
  itk::MetaDataDictionary       inputMetaDataDictionary = resampleImage->GetMetaDataDictionary();
  GenericTransformType::Pointer baseTransform = NULL;
  if( inputTransform == "ID" )
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
              << " wasn't of the expected type itk::VersorRigid3DTransform<double"
              << std::endl;
    return EXIT_FAILURE;
    }

  // Get measurement frame and its inverse from DWI scan
  std::vector<std::vector<double> > msrFrame;
  itk::ExposeMetaData<std::vector<std::vector<double> > >( inputMetaDataDictionary, "NRRD_measurement frame",
                                                           msrFrame );

  vnl_matrix_fixed<double, 3, 3> DWIMeasurementFrame;
  for( unsigned int i = 0; i < 3; i++ )
    {
    for( unsigned int j = 0; j < 3; j++ )
      {
      DWIMeasurementFrame[i][j] = msrFrame[i][j];
      }
    }
  std::cout << "DWI measurement frame (DWIMeasurementFrame): " << DWIMeasurementFrame << std::endl;
  vnl_matrix_fixed<double, 3, 3> DWIInverseMeasurementFrame = vnl_inverse( DWIMeasurementFrame );
  std::cout << "DWI inverse measurement frame (DWIInverseMeasurementFrame): " << DWIInverseMeasurementFrame
            << std::endl;

  // Resample DWI in place
  resampleImage = SetVectorImageRigidTransformInPlace<NrrdImageType>(rigidTransform.GetPointer(), resampleImage);

  std::cout << "Rigid transform matrix: " << rigidTransform->GetMatrix().GetVnlMatrix() << std::endl;
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
    std::cout << "Current gradient direction (before rotation): " << curGradientDirection << std::endl;
    RigidTransformType::Pointer inverseRigidTransform = RigidTransformType::New();
    const bool                  invertWorked = rigidTransform->GetInverse( inverseRigidTransform );
    curGradientDirection = inverseRigidTransform->GetMatrix().GetVnlMatrix() * DWIInverseMeasurementFrame
      * curGradientDirection;
    std::cout << "Current gradient direction (after rotation): " << curGradientDirection << std::endl;

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
    }

  // Set DWI measurement frame to identity by multiplying by its inverse
  // Update Image MetaData Dictionary with new measurement frame
  vnl_matrix_fixed<double, 3, 3> newMeasurementFrame = DWIInverseMeasurementFrame * DWIMeasurementFrame;
  std::cout << "New measurement frame (newMeasurementFrame): " << newMeasurementFrame << std::endl;
  std::vector<std::vector<double> > newMf(3);
  for( unsigned int i = 0; i < 3; i++ )
    {
    newMf[i].resize(3);
    for( unsigned int j = 0; j < 3; j++ )
      {
      newMf[i][j] = newMeasurementFrame[i][j];
      }
    }
  itk::EncapsulateMetaData<std::vector<std::vector<double> > >(
    resampleImage->GetMetaDataDictionary(), "NRRD_measurement frame", newMf );

  // Write out resampled in place DWI
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
  catch( itk::ExceptionObject & e )
    {
    std::cout << e << std::endl;
    }
  return EXIT_SUCCESS;
}
