#include <itkImage.h>
#include <itkImageFileWriter.h>
#include <itkImageFileReader.h>
#include <itkBinaryThresholdImageFilter.h>
#include <itkAndImageFilter.h>
#include <itkOrImageFilter.h>
#include <itkXorImageFilter.h>
#include <itkMedianImageFilter.h>
#include <itkMaximumImageFilter.h>
#include <itkBinaryErodeImageFilter.h>
#include <itkGradientAnisotropicDiffusionImageFilter.h>
#include <itkMaskImageFilter.h>
#include <itkBinaryDilateImageFilter.h>
#include <itkBinaryBallStructuringElement.h>


int main (int argc, char **argv)
{
  typedef itk::Image<unsigned char, 3>         ImageType;
  typedef itk::ImageFileReader<ImageType>      ImageReaderType;
  typedef itk::ImageFileWriter<ImageType>      ImageWriterType;

  
  ImageReaderType::Pointer classImageReader = ImageReaderType::New();
  classImageReader->SetFileName(argv[1]);
  classImageReader->Update();
  
  ImageReaderType::Pointer ventImageReader = ImageReaderType::New();
  ventImageReader->SetFileName(argv[2]);
  ventImageReader->Update();
  
  ImageReaderType::Pointer rCaudImageReader = ImageReaderType::New();
  rCaudImageReader->SetFileName(argv[3]);
  rCaudImageReader->Update();
  
  ImageReaderType::Pointer lCaudImageReader = ImageReaderType::New();
  lCaudImageReader->SetFileName(argv[4]);
  lCaudImageReader->Update();
  
  ImageReaderType::Pointer rPutImageReader = ImageReaderType::New();
  rPutImageReader->SetFileName(argv[5]);
  rPutImageReader->Update();
  
  ImageReaderType::Pointer lPutImageReader = ImageReaderType::New();
  lPutImageReader->SetFileName(argv[6]);
  lPutImageReader->Update();
  
  ImageReaderType::Pointer rThalImageReader = ImageReaderType::New();
  rThalImageReader->SetFileName(argv[7]);
  rThalImageReader->Update();
  
  ImageReaderType::Pointer lThalImageReader = ImageReaderType::New();
  lThalImageReader->SetFileName(argv[8]);
  lThalImageReader->Update();
  
  ImageReaderType::Pointer brainImageReader = ImageReaderType::New();
  brainImageReader->SetFileName(argv[9]);
  brainImageReader->Update();
  
  
  
  typedef itk::BinaryThresholdImageFilter<ImageType, ImageType>      BinaryFilterType;
  BinaryFilterType::Pointer binaryClassFilter = BinaryFilterType::New();
  binaryClassFilter->SetInput( classImageReader->GetOutput() );
  binaryClassFilter->SetLowerThreshold(9);
  binaryClassFilter->SetUpperThreshold(70);
  binaryClassFilter->SetOutsideValue(0);
  binaryClassFilter->SetInsideValue(1);
  binaryClassFilter->Update( );
  
  BinaryFilterType::Pointer binaryVentFilter = BinaryFilterType::New();
  binaryVentFilter->SetInput( ventImageReader->GetOutput() );
  binaryVentFilter->SetLowerThreshold(1);
  binaryVentFilter->SetUpperThreshold(255);
  binaryVentFilter->SetOutsideValue(0);
  binaryVentFilter->SetInsideValue(1);
  binaryVentFilter->Update( );
  
  BinaryFilterType::Pointer binaryRightCaudateFilter = BinaryFilterType::New();
  binaryRightCaudateFilter->SetInput( rCaudImageReader->GetOutput() );
  binaryRightCaudateFilter->SetLowerThreshold(1);
  binaryRightCaudateFilter->SetUpperThreshold(255);
  binaryRightCaudateFilter->SetOutsideValue(0);
  binaryRightCaudateFilter->SetInsideValue(1);
  binaryRightCaudateFilter->Update( );
  
  BinaryFilterType::Pointer binaryLeftCaudateFilter = BinaryFilterType::New();
  binaryLeftCaudateFilter->SetInput( lCaudImageReader->GetOutput() );
  binaryLeftCaudateFilter->SetLowerThreshold(1);
  binaryLeftCaudateFilter->SetUpperThreshold(255);
  binaryLeftCaudateFilter->SetOutsideValue(0);
  binaryLeftCaudateFilter->SetInsideValue(1);
  binaryLeftCaudateFilter->Update( );
  
  
  typedef itk::OrImageFilter<ImageType, ImageType, ImageType>   OrFilterType;
  OrFilterType::Pointer caudateFilter = OrFilterType::New();
  caudateFilter->SetInput1(binaryRightCaudateFilter->GetOutput());
  caudateFilter->SetInput2(binaryLeftCaudateFilter->GetOutput());
  caudateFilter->Update( );
  
  BinaryFilterType::Pointer binaryRightPutamenFilter = BinaryFilterType::New();
  binaryRightPutamenFilter->SetInput( rPutImageReader->GetOutput() );
  binaryRightPutamenFilter->SetLowerThreshold(1);
  binaryRightPutamenFilter->SetUpperThreshold(255);
  binaryRightPutamenFilter->SetOutsideValue(0);
  binaryRightPutamenFilter->SetInsideValue(1);
  binaryRightPutamenFilter->Update( );
  
  BinaryFilterType::Pointer binaryLeftPutamenFilter = BinaryFilterType::New();
  binaryLeftPutamenFilter->SetInput( lPutImageReader->GetOutput() );
  binaryLeftPutamenFilter->SetLowerThreshold(1);
  binaryLeftPutamenFilter->SetUpperThreshold(255);
  binaryLeftPutamenFilter->SetOutsideValue(0);
  binaryLeftPutamenFilter->SetInsideValue(1);
  binaryLeftPutamenFilter->Update( );
  
  OrFilterType::Pointer putamenFilter = OrFilterType::New();
  putamenFilter->SetInput1(binaryRightPutamenFilter->GetOutput());
  putamenFilter->SetInput2(binaryLeftPutamenFilter->GetOutput());
  putamenFilter->Update( );
  
  OrFilterType::Pointer basalGangliaFilter = OrFilterType::New();
  basalGangliaFilter->SetInput1(caudateFilter->GetOutput());
  basalGangliaFilter->SetInput2(putamenFilter->GetOutput());
  basalGangliaFilter->Update( );
  
  
  BinaryFilterType::Pointer binaryRightThalamusFilter = BinaryFilterType::New();
  binaryRightThalamusFilter->SetInput( rThalImageReader->GetOutput() );
  binaryRightThalamusFilter->SetLowerThreshold(1);
  binaryRightThalamusFilter->SetUpperThreshold(255);
  binaryRightThalamusFilter->SetOutsideValue(0);
  binaryRightThalamusFilter->SetInsideValue(1);
  binaryRightThalamusFilter->Update( );
  
  BinaryFilterType::Pointer binaryLeftThalamusFilter = BinaryFilterType::New();
  binaryLeftThalamusFilter->SetInput( lThalImageReader->GetOutput() );
  binaryLeftThalamusFilter->SetLowerThreshold(1);
  binaryLeftThalamusFilter->SetUpperThreshold(255);
  binaryLeftThalamusFilter->SetOutsideValue(0);
  binaryLeftThalamusFilter->SetInsideValue(1);
  binaryLeftThalamusFilter->Update( );
  
  OrFilterType::Pointer thalamusFilter = OrFilterType::New();
  thalamusFilter->SetInput1(binaryRightThalamusFilter->GetOutput());
  thalamusFilter->SetInput2(binaryLeftThalamusFilter->GetOutput());
  thalamusFilter->Update( );
  
  
  OrFilterType::Pointer subcorticalFilter = OrFilterType::New();
  subcorticalFilter->SetInput1(basalGangliaFilter->GetOutput());
  subcorticalFilter->SetInput2(thalamusFilter->GetOutput());
  subcorticalFilter->Update( );
  
  typedef itk::AndImageFilter<ImageType, ImageType, ImageType>   AndFilterType;
  AndFilterType::Pointer andFilter = AndFilterType::New();
  andFilter->SetInput1( binaryVentFilter->GetOutput() );
  andFilter->SetInput2( binaryClassFilter->GetOutput() );
  andFilter->Update();
  
  typedef itk::BinaryBallStructuringElement< unsigned char, 3 > KernelType;
  KernelType ballKernel;
  ballKernel.SetRadius(1);
  
  //unsigned long radius[3] = {1,1,1};
  
  typedef itk::BinaryDilateImageFilter< ImageType, ImageType, KernelType> DilateFilterType;
  DilateFilterType::Pointer dilateFilter = DilateFilterType::New();
  dilateFilter->SetInput( andFilter->GetOutput() );
  dilateFilter->SetKernel(ballKernel);
  dilateFilter->SetForegroundValue(1);
  dilateFilter->SetBackgroundValue(0);
  dilateFilter->SetDilateValue(1);
  dilateFilter->Update();
  
  OrFilterType::Pointer fillFilter = OrFilterType::New();
  fillFilter->SetInput1(dilateFilter->GetOutput());
  fillFilter->SetInput2(subcorticalFilter->GetOutput());
  fillFilter->Update( );
  
  ImageWriterType::Pointer writeResultFilter = ImageWriterType::New();
  writeResultFilter->SetInput( fillFilter->GetOutput() );
  writeResultFilter->SetFileName( argv[10] );
  writeResultFilter->Update();

  BinaryFilterType::Pointer rescaleVentFilter = BinaryFilterType::New();
  rescaleVentFilter->SetInput( fillFilter->GetOutput() );
  rescaleVentFilter->SetLowerThreshold(1);
  rescaleVentFilter->SetUpperThreshold(255);
  rescaleVentFilter->SetOutsideValue(0);
  rescaleVentFilter->SetInsideValue(230);
  rescaleVentFilter->Update( );
  
  typedef itk::MaximumImageFilter< ImageType, ImageType, ImageType > MaximumFilterType;
  MaximumFilterType::Pointer maximumFilter =  MaximumFilterType::New();
  maximumFilter->SetInput1( rescaleVentFilter->GetOutput() );
  maximumFilter->SetInput2( classImageReader->GetOutput() );
  maximumFilter->Update();
  
  typedef itk::MedianImageFilter<ImageType, ImageType> MedianFilterType;
  itk::Size<3> radius; 
  radius.Fill(1);
  MedianFilterType::Pointer medianFilter = MedianFilterType::New();
  medianFilter->SetInput( maximumFilter->GetOutput() );
  medianFilter->SetRadius(radius);
  medianFilter->Update();
  
  typedef itk::GradientAnisotropicDiffusionImageFilter<ImageType, ImageType> AnisotropicDiffusionType;
  AnisotropicDiffusionType::Pointer anisoFilter = AnisotropicDiffusionType::New();
  anisoFilter->SetInput( medianFilter->GetOutput() );
  anisoFilter->SetNumberOfIterations(5);
  anisoFilter->SetTimeStep(0.1);
  anisoFilter->SetConductanceParameter(1.0);
  anisoFilter->Update();
  
  ImageWriterType::Pointer writeClassFilter = ImageWriterType::New();
  writeClassFilter->SetInput( anisoFilter->GetOutput() );
  writeClassFilter->SetFileName( argv[11] );
  writeClassFilter->Update();

  /* Fix the Brain Mask */
  BinaryFilterType::Pointer binaryBrainFilter = BinaryFilterType::New();
  binaryBrainFilter->SetInput( brainImageReader->GetOutput() );
  binaryBrainFilter->SetLowerThreshold(1);
  binaryBrainFilter->SetUpperThreshold(255);
  binaryBrainFilter->SetOutsideValue(0);
  binaryBrainFilter->SetInsideValue(1);
  binaryBrainFilter->Update( );
  
  BinaryFilterType::Pointer binaryWhiteFilter = BinaryFilterType::New();
  binaryWhiteFilter->SetInput( anisoFilter->GetOutput() );
  binaryWhiteFilter->SetLowerThreshold(190);
  binaryWhiteFilter->SetUpperThreshold(255);
  binaryWhiteFilter->SetOutsideValue(0);
  binaryWhiteFilter->SetInsideValue(1);
  binaryWhiteFilter->Update( );
  
  typedef itk::BinaryErodeImageFilter< ImageType, ImageType, KernelType> ErodeFilterType;
  ErodeFilterType::Pointer erodeFilter = ErodeFilterType::New();
  erodeFilter->SetInput( binaryBrainFilter->GetOutput() );
  erodeFilter->SetKernel(ballKernel);
  erodeFilter->SetForegroundValue(1);
  erodeFilter->SetBackgroundValue(0);
  erodeFilter->SetErodeValue(1);
  erodeFilter->Update();
  
  typedef itk::XorImageFilter< ImageType, ImageType, ImageType> XorFilterType;
  XorFilterType::Pointer xorBrainFilter = XorFilterType::New();
  xorBrainFilter->SetInput1( binaryBrainFilter->GetOutput() );
  xorBrainFilter->SetInput2( erodeFilter->GetOutput() );
  xorBrainFilter->Update();
  
  XorFilterType::Pointer clipBrainFilter = XorFilterType::New();
  clipBrainFilter->SetInput1( xorBrainFilter->GetOutput() );
  clipBrainFilter->SetInput2( binaryBrainFilter->GetOutput() );
  clipBrainFilter->Update();
  
  ImageWriterType::Pointer writeBrainFilter = ImageWriterType::New();
  writeBrainFilter->SetInput( clipBrainFilter->GetOutput() );
  writeBrainFilter->SetFileName( argv[12] );
  writeBrainFilter->Update();
  
}

