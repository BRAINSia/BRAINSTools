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
#include <iostream>
#include <BRAINSFitHelper.h>
#include <itkEllipseSpatialObject.h>
#include <itkSpatialObjectToImageFilter.h>
#include <itkVersorRigid3DTransform.h>
#include <itkImageMaskSpatialObject.h>
#include <itkCastImageFilter.h>
#include <itkIO.h>

int main(int, char * *)
{
  using MaskPixelType = unsigned char;
  using PixelType = float;
  using ImageType = itk::Image<PixelType, 3>;
  using LocalMaskImageType = itk::Image<MaskPixelType, 3>;
  using EllipseSOType = itk::EllipseSpatialObject<3>;
  using SOToImageFilter = itk::SpatialObjectToImageFilter<EllipseSOType, ImageType>;
  using HelperType = itk::BRAINSFitHelper;
  using TransformType = EllipseSOType::TransformType;
  using ImageMaskSOType = itk::ImageMaskSpatialObject<3>;
  using CastType = itk::CastImageFilter<ImageType, LocalMaskImageType>;
  using CompositeTransformType = itk::CompositeTransform<double, 3>;

  // create two empty images
  ImageType::Pointer image1 = ImageType::New(),
    image2 = ImageType::New();

  ImageType::RegionType region;
  ImageType::SizeType   size;
  ImageType::IndexType  origin;
  size.Fill(100);
  origin.Fill(0);
  region.SetSize(size);
  region.SetIndex(origin);

  image1->SetRegions(region);
  image2->SetRegions(region);
  image1->Allocate();
  image2->Allocate();

  image1->FillBuffer(0);
  image2->FillBuffer(0);

  EllipseSOType::Pointer ellipse = EllipseSOType::New();
  TransformType::Pointer tfm = TransformType::New();
  tfm->SetIdentity();

  TransformType::OutputVectorType rotAxis;
  TransformType::OutputVectorType transVector;

  // and two ellipses, one of which is rotated and translated
  EllipseSOType::ArrayType ePar;
  ePar[0] = 10;
  ePar[1] = 20;
  ePar[2] = 40;

  transVector.Fill(50);

  tfm->Translate(transVector);
  ellipse->SetRadius(ePar);
  ellipse->SetObjectToWorldTransform(tfm);

  // convert ellipses to binary images
  SOToImageFilter::Pointer e2image = SOToImageFilter::New();
  SOToImageFilter::Pointer etfm2image = SOToImageFilter::New();
  e2image->SetInput(ellipse);
  e2image->SetSize(size);
  e2image->Update();
  ImageType::Pointer eImage = e2image->GetOutput();

  rotAxis.Fill(1.);
  float rotAngle = 3.14 / 3.;
  tfm->Rotate3D(rotAxis, rotAngle);
  transVector[0] = 10;
  transVector[1] = -5;
  transVector[2] = 15;
  tfm->Translate(transVector);

  // translate the second ellipse by a different tfm and rotate
  ellipse->SetObjectToWorldTransform(tfm);
  ellipse->Initialize();

  etfm2image->SetInput(ellipse);
  etfm2image->SetSize(size);
  etfm2image->Update();

  ImageType::Pointer eTfmImage = etfm2image->GetOutput();

  using VersorRigid3DTransformType = itk::VersorRigid3DTransform<double>;
  VersorRigid3DTransformType::Pointer tempCopy = VersorRigid3DTransformType::New();

  // images and masks passed to helper are identical, but only masks will be used
  CastType::Pointer fixedImageCast = CastType::New();
  CastType::Pointer movingImageCast = CastType::New();

  ImageMaskSOType::Pointer fixedMask = ImageMaskSOType::New();
  ImageMaskSOType::Pointer movingMask = ImageMaskSOType::New();

  fixedImageCast->SetInput(eImage);
  movingImageCast->SetInput(eTfmImage);

  // need to create spatial objects back from binary images
  fixedMask->SetImage(fixedImageCast->GetOutput() );
  fixedMask->Update(); // Replaced old ComputeObjectToWorldTransform with new Update()
  movingMask->SetImage(movingImageCast->GetOutput() );
  movingMask->Update(); // Replaced old ComputeObjectToWorldTransform with new Update()

  std::vector<std::string> transformTypeVector;

  // initialize the helper and run
  HelperType::Pointer myHelper = HelperType::New();
  myHelper->SetFixedVolume(eImage);
  myHelper->SetMovingVolume(eTfmImage);
  myHelper->SetFixedBinaryVolume(fixedMask);
  myHelper->SetMovingBinaryVolume(movingMask);
  myHelper->SetCurrentGenericTransform(nullptr);
  myHelper->SetInitializeTransformMode("useCenterOfROIAlign");
  myHelper->SetTransformType(transformTypeVector);
  myHelper->Update();

  using GenericTransformType = itk::Transform<double, 3, 3>;
  GenericTransformType::Pointer currentGenericTransform = myHelper->GetCurrentGenericTransform().GetPointer();

  const CompositeTransformType::ConstPointer genericCompositeTransform =
    dynamic_cast<const CompositeTransformType *>( currentGenericTransform.GetPointer() );
  if( genericCompositeTransform.IsNull() )
    {
    itkGenericExceptionMacro(<<"Error in transform type conversion");
    }

  VersorRigid3DTransformType::ConstPointer versor3D =
    dynamic_cast<const VersorRigid3DTransformType *>(genericCompositeTransform->GetNthTransform(0).GetPointer() );
  if( versor3D.IsNull() )
    {
    itkGenericExceptionMacro(<<"Error in transform type conversion");
    }

  // check translation vector -- should match the difference in the ellipse translation vectors
  TransformType::OutputVectorType recoveredTransVector = versor3D->GetTranslation();

  if( transVector != recoveredTransVector )
    {
    return EXIT_FAILURE;
    }

  return EXIT_SUCCESS;
}
