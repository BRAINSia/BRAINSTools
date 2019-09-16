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

/*
static std::string myGetEnv(std::string p)
{
  const char *  const temp = std::getenv(p.c_str());
  if (temp != nullptr)
  {
    return std::string(temp);
  }
  return std::string("");
}
*/

int
main(int, char **)
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

  ImageType::SizeType size;
  size.Fill(100);

  TransformType::Pointer tfm = TransformType::New();
  tfm->SetIdentity();

  EllipseSOType::Pointer ellipse = EllipseSOType::New();
  {
    // and two ellipses, one of which is rotated and translated
    EllipseSOType::ArrayType ePar;
    ePar[0] = 10;
    ePar[1] = 20;
    ePar[2] = 40;
    ellipse->SetRadiusInObjectSpace(ePar);
  }

  TransformType::OutputVectorType transVector;
  ImageType::Pointer              eImage;
  {
    transVector.Fill(50);
    tfm->Translate(transVector);
    ellipse->SetObjectToWorldTransform(tfm);
    ellipse->Initialize();
    // convert ellipses to binary images
    SOToImageFilter::Pointer e2image = SOToImageFilter::New();
    e2image->SetInput(ellipse);
    e2image->SetSize(size);
    e2image->Update();
    eImage = e2image->GetOutput();
  }

  ImageType::Pointer eTfmImage;
  {
    TransformType::OutputVectorType rotAxis;
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

    SOToImageFilter::Pointer etfm2image = SOToImageFilter::New();
    etfm2image->SetInput(ellipse);
    etfm2image->SetSize(size);
    etfm2image->Update();
    eTfmImage = etfm2image->GetOutput();
  }

  std::cout << eImage->GetSpacing() << std::endl;
  std::cout << eTfmImage->GetSpacing() << std::endl;

  using VersorRigid3DTransformType = itk::VersorRigid3DTransform<double>;
  VersorRigid3DTransformType::Pointer tempCopy = VersorRigid3DTransformType::New();

  // images and masks passed to helper are identical, but only masks will be used
  CastType::Pointer fixedImageCast = CastType::New();
  fixedImageCast->SetInput(eImage);
  CastType::Pointer movingImageCast = CastType::New();
  movingImageCast->SetInput(eTfmImage);

  // std::string prefix=myGetEnv("DEBUG_PREFIX")
  itkUtil::WriteImage<LocalMaskImageType>(fixedImageCast->GetOutput(), "/tmp/fixed.nii.gz");
  itkUtil::WriteImage<LocalMaskImageType>(movingImageCast->GetOutput(), "/tmp/moving.nii.gz");

  // need to create spatial objects back from binary images
  ImageMaskSOType::Pointer fixedMask = ImageMaskSOType::New();
  fixedMask->SetImage(fixedImageCast->GetOutput());
  fixedMask->Update(); // Replaced old ComputeObjectToWorldTransform with new Update()

  ImageMaskSOType::Pointer movingMask = ImageMaskSOType::New();
  movingMask->SetImage(movingImageCast->GetOutput());
  movingMask->Update(); // Replaced old ComputeObjectToWorldTransform with new Update()

  std::vector<std::string> transformTypeVector;
  transformTypeVector.push_back(std::string("Rigid"));

  // initialize the helper and run
  HelperType::Pointer myHelper = HelperType::New();
  myHelper->SetFixedVolume(eImage);
  myHelper->SetCostMetricName("MSE"); // MSE, images are binary, and MMI does not work on binary images!
  myHelper->SetMovingVolume(eTfmImage);
  myHelper->SetFixedBinaryVolume(fixedMask);
  myHelper->SetMovingBinaryVolume(movingMask);
  myHelper->SetCurrentGenericTransform(nullptr);
  myHelper->SetInitializeTransformMode("useCenterOfROIAlign");
  myHelper->SetTransformType(transformTypeVector);
  // std::vector<double> minStepLength;
  // minStepLength.push_back(0.001);
  // myHelper->SetMinimumStepLength(minStepLength);
  myHelper->SetDebugLevel(100);
  myHelper->PrintCommandLine(true, "DEBUGTESTFAILURES");
  myHelper->Update();


  using GenericTransformType = itk::Transform<double, 3, 3>;
  GenericTransformType::Pointer currentGenericTransform = myHelper->GetCurrentGenericTransform().GetPointer();

  const CompositeTransformType::ConstPointer genericCompositeTransform =
    dynamic_cast<const CompositeTransformType *>(currentGenericTransform.GetPointer());
  if (genericCompositeTransform.IsNull())
  {
    itkGenericExceptionMacro(<< "Error in transform type conversion");
  }

  VersorRigid3DTransformType::ConstPointer versor3D =
    dynamic_cast<const VersorRigid3DTransformType *>(genericCompositeTransform->GetNthTransform(0).GetPointer());
  if (versor3D.IsNull())
  {
    itkGenericExceptionMacro(<< "Error in transform type conversion");
  }

  // check translation vector -- should match the difference in the ellipse translation vectors
  TransformType::OutputVectorType recoveredTransVector = versor3D->GetTranslation();

  double error = 0.0;
  for (size_t i = 0; i < 3; ++i)
  {
    const double e = (transVector[i] - recoveredTransVector[i]);
    error += e * e;
  }
  constexpr double tolerance = 0.05;
  std::cout << "PLOT," << error << "," << transVector << "," << recoveredTransVector << "," << tolerance << std::endl;
  if (error > tolerance)
  {
    std::cout << "ERROR term too big: " << error << " must be less than " << tolerance << std::endl;
    std::cout << "Known Translation Vector: " << transVector << std::endl;
    std::cout << "Recovered Vector        : " << recoveredTransVector << std::endl;
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
