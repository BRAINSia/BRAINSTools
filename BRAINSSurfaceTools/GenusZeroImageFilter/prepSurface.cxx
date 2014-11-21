#include <itkImage.h>
#include <itkImageFileWriter.h>
#include <itkImageFileReader.h>
#include <itkBinaryThresholdImageFilter.h>
#include <itkAndImageFilter.h>
#include <itkOrImageFilter.h>
#include <itkNotImageFilter.h>
#include <itkXorImageFilter.h>
#include <itkMedianImageFilter.h>
#include <itkMaximumImageFilter.h>
#include <itkBinaryErodeImageFilter.h>
#include <itkGradientAnisotropicDiffusionImageFilter.h>
#include <itkMaskImageFilter.h>
#include <itkCastImageFilter.h>
#include <itkBinaryDilateImageFilter.h>
#include <itkBinaryBallStructuringElement.h>
#include <itkRelabelComponentImageFilter.h>
#include <itkConnectedComponentImageFilter.h>

#include "prepSurfaceCLP.h"

int main(int argc, char * *argv)
{
  PARSE_ARGS;

  typedef itk::Image<unsigned char, 3>    ImageType;
  typedef itk::Image<float, 3>            FloatImageType;
  typedef itk::ImageFileReader<ImageType> ImageReaderType;
  typedef itk::ImageFileWriter<ImageType> ImageWriterType;

  ImageReaderType::Pointer classImageReader = ImageReaderType::New();
  classImageReader->SetFileName( tissueClassVolume );
  classImageReader->Update();

  ImageReaderType::Pointer ventImageReader = ImageReaderType::New();
  ventImageReader->SetFileName(ventricleVolume);
  ventImageReader->Update();

  ImageReaderType::Pointer rCaudImageReader = ImageReaderType::New();
  rCaudImageReader->SetFileName(rightCaudateVolume);
  rCaudImageReader->Update();

  ImageReaderType::Pointer lCaudImageReader = ImageReaderType::New();
  lCaudImageReader->SetFileName(leftCaudateVolume);
  lCaudImageReader->Update();

  ImageReaderType::Pointer rPutImageReader = ImageReaderType::New();
  rPutImageReader->SetFileName(rightPutamenVolume);
  rPutImageReader->Update();

  ImageReaderType::Pointer lPutImageReader = ImageReaderType::New();
  lPutImageReader->SetFileName(leftPutamenVolume);
  lPutImageReader->Update();

  ImageReaderType::Pointer rThalImageReader = ImageReaderType::New();
  rThalImageReader->SetFileName(rightThalamusVolume);
  rThalImageReader->Update();

  ImageReaderType::Pointer lThalImageReader = ImageReaderType::New();
  lThalImageReader->SetFileName(leftThalamusVolume);
  lThalImageReader->Update();

  ImageReaderType::Pointer brainImageReader = ImageReaderType::New();
  brainImageReader->SetFileName(brainVolume);
  brainImageReader->Update();

  ImageReaderType::Pointer leftHemisphereReader = ImageReaderType::New();
  leftHemisphereReader->SetFileName(leftHemisphereVolume);
  leftHemisphereReader->Update();

  ImageReaderType::Pointer rightHemisphereReader = ImageReaderType::New();
  rightHemisphereReader->SetFileName(rightHemisphereVolume);
  rightHemisphereReader->Update();

  ImageReaderType::Pointer clipImageReader = ImageReaderType::New();
  clipImageReader->SetFileName( clipVolume );
  clipImageReader->Update();

  typedef itk::BinaryThresholdImageFilter<ImageType, ImageType> BinaryFilterType;
  BinaryFilterType::Pointer binaryVentFilter = BinaryFilterType::New();
  binaryVentFilter->SetInput( ventImageReader->GetOutput() );
  binaryVentFilter->SetLowerThreshold(1);
  binaryVentFilter->SetUpperThreshold(255);
  binaryVentFilter->SetOutsideValue(0);
  binaryVentFilter->SetInsideValue(1);
  binaryVentFilter->Update();

  BinaryFilterType::Pointer binaryBrainFilter = BinaryFilterType::New();
  binaryBrainFilter->SetInput( brainImageReader->GetOutput() );
  binaryBrainFilter->SetLowerThreshold(1);
  binaryBrainFilter->SetUpperThreshold(255);
  binaryBrainFilter->SetOutsideValue(0);
  binaryBrainFilter->SetInsideValue(1);
  binaryBrainFilter->Update();

  BinaryFilterType::Pointer binaryClipFilter = BinaryFilterType::New();
  binaryClipFilter->SetInput( clipImageReader->GetOutput() );
  binaryClipFilter->SetLowerThreshold(1);
  binaryClipFilter->SetUpperThreshold(255);
  binaryClipFilter->SetOutsideValue(0);
  binaryClipFilter->SetInsideValue(1);
  binaryClipFilter->Update();

  // Invert clip region
  typedef itk::NotImageFilter<ImageType, ImageType> NotFilterType;
  NotFilterType::Pointer invertFilter = NotFilterType::New();
  invertFilter->SetInput( binaryClipFilter->GetOutput() );
  invertFilter->Update();

  // Combine Brain region and inverse of clip region
  typedef itk::AndImageFilter<ImageType, ImageType, ImageType> AndFilterType;
  AndFilterType::Pointer brainRegionFilter = AndFilterType::New();
  brainRegionFilter->SetInput1(binaryBrainFilter->GetOutput() );
  brainRegionFilter->SetInput2(invertFilter->GetOutput() );
  brainRegionFilter->Update();

  // Generate Caudate Region
  BinaryFilterType::Pointer binaryRightCaudateFilter = BinaryFilterType::New();
  binaryRightCaudateFilter->SetInput( rCaudImageReader->GetOutput() );
  binaryRightCaudateFilter->SetLowerThreshold(1);
  binaryRightCaudateFilter->SetUpperThreshold(255);
  binaryRightCaudateFilter->SetOutsideValue(0);
  binaryRightCaudateFilter->SetInsideValue(1);
  binaryRightCaudateFilter->Update();

  BinaryFilterType::Pointer binaryLeftCaudateFilter = BinaryFilterType::New();
  binaryLeftCaudateFilter->SetInput( lCaudImageReader->GetOutput() );
  binaryLeftCaudateFilter->SetLowerThreshold(1);
  binaryLeftCaudateFilter->SetUpperThreshold(255);
  binaryLeftCaudateFilter->SetOutsideValue(0);
  binaryLeftCaudateFilter->SetInsideValue(1);
  binaryLeftCaudateFilter->Update();

  typedef itk::OrImageFilter<ImageType, ImageType, ImageType> OrFilterType;
  OrFilterType::Pointer caudateFilter = OrFilterType::New();
  caudateFilter->SetInput1(binaryRightCaudateFilter->GetOutput() );
  caudateFilter->SetInput2(binaryLeftCaudateFilter->GetOutput() );
  caudateFilter->Update();

  ImageWriterType::Pointer writeTmpFilter1 = ImageWriterType::New();
  writeTmpFilter1->SetInput( caudateFilter->GetOutput() );
  writeTmpFilter1->SetFileName( "caudate.nii.gz" );
  writeTmpFilter1->Update();

  // Generate Putamen Region
  BinaryFilterType::Pointer binaryRightPutamenFilter = BinaryFilterType::New();
  binaryRightPutamenFilter->SetInput( rPutImageReader->GetOutput() );
  binaryRightPutamenFilter->SetLowerThreshold(1);
  binaryRightPutamenFilter->SetUpperThreshold(255);
  binaryRightPutamenFilter->SetOutsideValue(0);
  binaryRightPutamenFilter->SetInsideValue(1);
  binaryRightPutamenFilter->Update();

  BinaryFilterType::Pointer binaryLeftPutamenFilter = BinaryFilterType::New();
  binaryLeftPutamenFilter->SetInput( lPutImageReader->GetOutput() );
  binaryLeftPutamenFilter->SetLowerThreshold(1);
  binaryLeftPutamenFilter->SetUpperThreshold(255);
  binaryLeftPutamenFilter->SetOutsideValue(0);
  binaryLeftPutamenFilter->SetInsideValue(1);
  binaryLeftPutamenFilter->Update();

  OrFilterType::Pointer putamenFilter = OrFilterType::New();
  putamenFilter->SetInput1(binaryRightPutamenFilter->GetOutput() );
  putamenFilter->SetInput2(binaryLeftPutamenFilter->GetOutput() );
  putamenFilter->Update();

  ImageWriterType::Pointer writeTmpFilter2 = ImageWriterType::New();
  writeTmpFilter2->SetInput( putamenFilter->GetOutput() );
  writeTmpFilter2->SetFileName( "putamen.nii.gz" );
  writeTmpFilter2->Update();

  // Basal Ganglia =  Caudate + Putamen
  OrFilterType::Pointer basalGangliaFilter = OrFilterType::New();
  basalGangliaFilter->SetInput1(caudateFilter->GetOutput() );
  basalGangliaFilter->SetInput2(putamenFilter->GetOutput() );
  basalGangliaFilter->Update();

  ImageWriterType::Pointer writeTmpFilter3 = ImageWriterType::New();
  writeTmpFilter3->SetInput( basalGangliaFilter->GetOutput() );
  writeTmpFilter3->SetFileName( "bg.nii.gz" );
  writeTmpFilter3->Update();

  // Generate Thalamus Region
  BinaryFilterType::Pointer binaryRightThalamusFilter = BinaryFilterType::New();
  binaryRightThalamusFilter->SetInput( rThalImageReader->GetOutput() );
  binaryRightThalamusFilter->SetLowerThreshold(1);
  binaryRightThalamusFilter->SetUpperThreshold(255);
  binaryRightThalamusFilter->SetOutsideValue(0);
  binaryRightThalamusFilter->SetInsideValue(1);
  binaryRightThalamusFilter->Update();

  BinaryFilterType::Pointer binaryLeftThalamusFilter = BinaryFilterType::New();
  binaryLeftThalamusFilter->SetInput( lThalImageReader->GetOutput() );
  binaryLeftThalamusFilter->SetLowerThreshold(1);
  binaryLeftThalamusFilter->SetUpperThreshold(255);
  binaryLeftThalamusFilter->SetOutsideValue(0);
  binaryLeftThalamusFilter->SetInsideValue(1);
  binaryLeftThalamusFilter->Update();

  OrFilterType::Pointer thalamusFilter = OrFilterType::New();
  thalamusFilter->SetInput1(binaryRightThalamusFilter->GetOutput() );
  thalamusFilter->SetInput2(binaryLeftThalamusFilter->GetOutput() );
  thalamusFilter->Update();

  ImageWriterType::Pointer writeTmpFilter4 = ImageWriterType::New();
  writeTmpFilter4->SetInput( thalamusFilter->GetOutput() );
  writeTmpFilter4->SetFileName( "thalamus.nii.gz" );
  writeTmpFilter4->Update();

  // Subcortical = Basal Ganglia + Thalamus
  OrFilterType::Pointer subcorticalFilter = OrFilterType::New();
  subcorticalFilter->SetInput1(basalGangliaFilter->GetOutput() );
  subcorticalFilter->SetInput2(thalamusFilter->GetOutput() );
  subcorticalFilter->Update();

  ImageWriterType::Pointer writeTmpFilter5 = ImageWriterType::New();
  writeTmpFilter5->SetInput( subcorticalFilter->GetOutput() );
  writeTmpFilter5->SetFileName( "subcortical.nii.gz" );
  writeTmpFilter5->Update();

  // Subcortical + Ventricles
  OrFilterType::Pointer fillFilter = OrFilterType::New();
  fillFilter->SetInput1(binaryVentFilter->GetOutput() );
  fillFilter->SetInput2(subcorticalFilter->GetOutput() );
  fillFilter->Update();

  ImageWriterType::Pointer writeTmpFilter6 = ImageWriterType::New();
  writeTmpFilter6->SetInput( fillFilter->GetOutput() );
  writeTmpFilter6->SetFileName( "fill.nii.gz" );
  writeTmpFilter6->Update();

  BinaryFilterType::Pointer rescaleVentFilter = BinaryFilterType::New();
  rescaleVentFilter->SetInput( fillFilter->GetOutput() );
  rescaleVentFilter->SetLowerThreshold(1);
  rescaleVentFilter->SetUpperThreshold(255);
  rescaleVentFilter->SetOutsideValue(0);
  rescaleVentFilter->SetInsideValue(230);
  rescaleVentFilter->Update();

  typedef itk::MaximumImageFilter<ImageType, ImageType, ImageType> MaximumFilterType;
  MaximumFilterType::Pointer maximumFilter =  MaximumFilterType::New();
  maximumFilter->SetInput1( rescaleVentFilter->GetOutput() );
  maximumFilter->SetInput2( classImageReader->GetOutput() );
  maximumFilter->Update();

  /* Now clip the image to the Hemispheres and Brain*/
  typedef itk::MaskImageFilter<ImageType, ImageType> MaskFilterType;
  MaskFilterType::Pointer maskLeftFilter = MaskFilterType::New();
  maskLeftFilter->SetInput1( maximumFilter->GetOutput() );
  maskLeftFilter->SetInput2( leftHemisphereReader->GetOutput() );
  maskLeftFilter->SetOutsideValue( 0 );
  maskLeftFilter->Update();

  MaskFilterType::Pointer maskRightFilter = MaskFilterType::New();
  maskRightFilter->SetInput1( maximumFilter->GetOutput() );
  maskRightFilter->SetInput2( rightHemisphereReader->GetOutput() );
  maskRightFilter->SetOutsideValue( 0 );
  maskLeftFilter->Update();

  MaskFilterType::Pointer maskLeftClipFilter = MaskFilterType::New();
  maskLeftClipFilter->SetInput1( maskLeftFilter->GetOutput() );
  maskLeftClipFilter->SetInput2( brainRegionFilter->GetOutput() );
  maskLeftClipFilter->SetOutsideValue( 0 );
  maskLeftClipFilter->Update();

  MaskFilterType::Pointer maskRightClipFilter = MaskFilterType::New();
  maskRightClipFilter->SetInput1( maskRightFilter->GetOutput() );
  maskRightClipFilter->SetInput2( brainRegionFilter->GetOutput() );
  maskRightClipFilter->SetOutsideValue( 0 );
  maskRightClipFilter->Update();

  /* Now Filter the image */
  ImageType::Pointer rightTissueClass = maskRightClipFilter->GetOutput();
  ImageType::Pointer leftTissueClass = maskLeftClipFilter->GetOutput();

  if( medianFilter )
    {
    typedef itk::MedianImageFilter<ImageType, ImageType> MedianFilterType;
    itk::Size<3> radius;
    radius.SetElement(0, medianFilterSize[0]);
    radius.SetElement(1, medianFilterSize[0]);
    radius.SetElement(2, medianFilterSize[0]);

    MedianFilterType::Pointer medianLeftFilter = MedianFilterType::New();
    medianLeftFilter->SetInput( maskLeftClipFilter->GetOutput() );
    medianLeftFilter->SetRadius(radius);
    medianLeftFilter->Update();
    leftTissueClass = medianLeftFilter->GetOutput();

    MedianFilterType::Pointer medianRightFilter = MedianFilterType::New();
    medianRightFilter->SetInput( maskRightClipFilter->GetOutput() );
    medianRightFilter->SetRadius(radius);
    medianRightFilter->Update();
    rightTissueClass = medianRightFilter->GetOutput();
    }

  if( anisoDiffusionFilter )
    {
    typedef itk::CastImageFilter<ImageType, FloatImageType> CastFloatType;
    CastFloatType::Pointer leftFloatCast = CastFloatType::New();
    leftFloatCast->SetInput( leftTissueClass );
    leftFloatCast->Update();

    typedef itk::GradientAnisotropicDiffusionImageFilter<FloatImageType, FloatImageType> AnisotropicDiffusionType;
    AnisotropicDiffusionType::Pointer anisoLeftFilter = AnisotropicDiffusionType::New();
    anisoLeftFilter->SetInput( leftFloatCast->GetOutput() );
    anisoLeftFilter->SetNumberOfIterations( iterations );
    anisoLeftFilter->SetTimeStep( static_cast<double>(stepSize) );
    anisoLeftFilter->SetConductanceParameter( conductance );
    anisoLeftFilter->Update();

    typedef itk::CastImageFilter<FloatImageType, ImageType> CastImageType;
    CastImageType::Pointer leftImageCast = CastImageType::New();
    leftImageCast->SetInput( anisoLeftFilter->GetOutput() );
    leftImageCast->Update();
    leftTissueClass = leftImageCast->GetOutput();

    CastFloatType::Pointer rightFloatCast = CastFloatType::New();
    rightFloatCast->SetInput( rightTissueClass );
    rightFloatCast->Update();

    AnisotropicDiffusionType::Pointer anisoRightFilter = AnisotropicDiffusionType::New();
    anisoRightFilter->SetInput( rightFloatCast->GetOutput() );
    anisoRightFilter->SetNumberOfIterations( iterations );
    anisoRightFilter->SetTimeStep( static_cast<double>(stepSize) );
    anisoRightFilter->SetConductanceParameter( conductance );
    anisoRightFilter->Update();

    CastImageType::Pointer rightImageCast = CastImageType::New();
    rightImageCast->SetInput( anisoRightFilter->GetOutput() );
    rightImageCast->Update();
    rightTissueClass = rightImageCast->GetOutput();
    }

  /* Create Binary Images for Surface Generation */
  BinaryFilterType::Pointer binaryLeftBrainFilter = BinaryFilterType::New();
  binaryLeftBrainFilter->SetInput( leftTissueClass );
  binaryLeftBrainFilter->SetLowerThreshold(190);
  binaryLeftBrainFilter->SetUpperThreshold(255);
  binaryLeftBrainFilter->SetOutsideValue(0);
  binaryLeftBrainFilter->SetInsideValue(1);
  binaryLeftBrainFilter->Update();

  typedef itk::ConnectedComponentImageFilter<ImageType, ImageType> ConnectedComponentFilterType;
  ConnectedComponentFilterType::Pointer leftConnectedRegions = ConnectedComponentFilterType::New();
  leftConnectedRegions->SetInput( binaryLeftBrainFilter->GetOutput() );
  leftConnectedRegions->Update();

  typedef itk::RelabelComponentImageFilter<ImageType, ImageType> RelabelImageFilterType;
  RelabelImageFilterType::Pointer relabelLeftFilter = RelabelImageFilterType::New();
  relabelLeftFilter->SetInput( leftConnectedRegions->GetOutput() );
  relabelLeftFilter->Update();

  BinaryFilterType::Pointer binaryLeftLabelFilter = BinaryFilterType::New();
  binaryLeftLabelFilter->SetInput( relabelLeftFilter->GetOutput() );
  binaryLeftLabelFilter->SetLowerThreshold(1);
  binaryLeftLabelFilter->SetUpperThreshold(1);
  binaryLeftLabelFilter->SetOutsideValue(0);
  binaryLeftLabelFilter->SetInsideValue(1);
  binaryLeftLabelFilter->Update();

  BinaryFilterType::Pointer binaryRightBrainFilter = BinaryFilterType::New();
  binaryRightBrainFilter->SetInput( rightTissueClass );
  binaryRightBrainFilter->SetLowerThreshold(190);
  binaryRightBrainFilter->SetUpperThreshold(255);
  binaryRightBrainFilter->SetOutsideValue(0);
  binaryRightBrainFilter->SetInsideValue(1);
  binaryRightBrainFilter->Update();

  ConnectedComponentFilterType::Pointer rightConnectedRegions = ConnectedComponentFilterType::New();
  rightConnectedRegions->SetInput( binaryRightBrainFilter->GetOutput() );
  rightConnectedRegions->Update();

  RelabelImageFilterType::Pointer relabelRightFilter = RelabelImageFilterType::New();
  relabelRightFilter->SetInput( rightConnectedRegions->GetOutput() );
  relabelRightFilter->Update();

  BinaryFilterType::Pointer binaryRightLabelFilter = BinaryFilterType::New();
  binaryRightLabelFilter->SetInput( relabelRightFilter->GetOutput() );
  binaryRightLabelFilter->SetLowerThreshold(1);
  binaryRightLabelFilter->SetUpperThreshold(1);
  binaryRightLabelFilter->SetOutsideValue(0);
  binaryRightLabelFilter->SetInsideValue(1);
  binaryRightLabelFilter->Update();

  /* Write the results */
  ImageWriterType::Pointer writeBinaryLeftFilter = ImageWriterType::New();
  writeBinaryLeftFilter->SetInput( binaryLeftLabelFilter->GetOutput() );
  writeBinaryLeftFilter->SetFileName( leftBinarySurfaceVolume );
  writeBinaryLeftFilter->Update();

  ImageWriterType::Pointer writeBinaryRightFilter = ImageWriterType::New();
  writeBinaryRightFilter->SetInput( binaryRightLabelFilter->GetOutput() );
  writeBinaryRightFilter->SetFileName( rightBinarySurfaceVolume );
  writeBinaryRightFilter->Update();

  ImageWriterType::Pointer writeTissueClassFilter = ImageWriterType::New();
  writeTissueClassFilter->SetInput( maximumFilter->GetOutput() );
  writeTissueClassFilter->SetFileName( filledClassVolume );
  writeTissueClassFilter->Update();

}
