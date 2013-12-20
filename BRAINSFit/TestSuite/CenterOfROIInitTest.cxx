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
  typedef unsigned char                                             MaskPixelType;
  typedef float                                                     PixelType;
  typedef itk::Image<PixelType, 3>                                  ImageType;
  typedef itk::Image<MaskPixelType, 3>                              LocalMaskImageType;
  typedef itk::EllipseSpatialObject<3>                              EllipseSOType;
  typedef itk::SpatialObjectToImageFilter<EllipseSOType, ImageType> SOToImageFilter;
  typedef itk::BRAINSFitHelper                                      HelperType;
  typedef EllipseSOType::TransformType                              TransformType;
  typedef itk::ImageMaskSpatialObject<3>                            ImageMaskSOType;
  typedef itk::CastImageFilter<ImageType, LocalMaskImageType>       CastType;

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
  fixedMask->ComputeObjectToWorldTransform();
  movingMask->SetImage(movingImageCast->GetOutput() );
  movingMask->ComputeObjectToWorldTransform();

  std::vector<std::string> transformTypeVector;

  // initialize the helper and run
  HelperType::Pointer myHelper = HelperType::New();
  myHelper->SetFixedVolume(eImage);
  myHelper->SetMovingVolume(eTfmImage);
  myHelper->SetFixedBinaryVolume(fixedMask);
  myHelper->SetMovingBinaryVolume(movingMask);
  myHelper->SetCurrentGenericTransform(NULL);
  myHelper->SetInitializeTransformMode("useCenterOfROIAlign");
  myHelper->SetTransformType(transformTypeVector);
  myHelper->Update();

  GenericTransformType::Pointer currentGenericTransform = myHelper->GetCurrentGenericTransform();

  VersorRigid3DTransformType::ConstPointer versor3D =
    dynamic_cast<const VersorRigid3DTransformType *>(currentGenericTransform.GetPointer() );

  // check translation vector -- should match the difference in the ellipse translation vectors
  TransformType::OutputVectorType recoveredTransVector = versor3D->GetTranslation();

  if( transVector != recoveredTransVector )
    {
    return EXIT_FAILURE;
    }

  return EXIT_SUCCESS;
}
