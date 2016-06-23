//
// Created by Leinoff, Alexander on 5/20/16.
//

#include <itkTransformToDisplacementFieldFilter.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <BRAINSRefacerCLP.h>
#include <Slicer3LandmarkIO.h>
#include <itkBinaryThresholdImageFilter.h>
#include <itkBSplineTransform.h>
#include <itkResampleImageFilter.h>
#include <itkDanielssonDistanceMapImageFilter.h>
#include <itkDivideImageFilter.h>
#include <itkVectorIndexSelectionCastImageFilter.h>
#include <itkTransformFileWriter.h>
#include <itkMultiplyImageFilter.h>
#include <itkBSplineTransformInitializer.h>
#include <itkComposeImageFilter.h>
#include <itkDisplacementFieldTransform.h>
#include <itkSubtractImageFilter.h>
#include <map>
#include <string>
#include <itkSignedMaurerDistanceMapImageFilter.h>
#include <itkThresholdImageFilter.h>

#include "CreateRandomBSpline.h"
#include "CombineBSplineWithDisplacement.h"
#include "MaskFromLandmarks.h"

//Convienience function to write images
template< typename TImageType >
void WriteImage(std::string filename, TImageType *image)
{
  std::cout << "Writing Image: " << filename << std::endl;
  typedef itk::ImageFileWriter<TImageType> FileWriterType;
  typename FileWriterType::Pointer fileWriter = FileWriterType::New();

  fileWriter->SetInput(image);
  fileWriter->SetFileName(filename);
  fileWriter->Update();
  std::cout << "\tdone writing Image: " << filename << std::endl;
}

//Convienience function to write images
template< typename TImageType >
void WriteSmartImage(std::string filename, typename TImageType::Pointer image)
{
  std::cout << "Writing Image: " << filename << std::endl;
  typedef itk::ImageFileWriter<TImageType> FileWriterType;
  typename FileWriterType::Pointer fileWriter = FileWriterType::New();

  fileWriter->SetInput(image);
  fileWriter->SetFileName(filename);
  fileWriter->Update();
  std::cout << "\tdone writing Image: " << filename << std::endl;
}


//Convienience function to write transforms
template< typename TTransformType >
void WriteTransform(std::string transformFileName, TTransformType transform )
{
  std::cout << "Writing Transform: " << transformFileName << std::endl;
  typedef itk::TransformFileWriter TransformWriterType;
  TransformWriterType::Pointer transformWriter = TransformWriterType::New();

  transformWriter->SetInput(transform);
  transformWriter->SetFileName(transformFileName);
  transformWriter->Update();
  std::cout << "\t done writing Transform: " << transformFileName << std::endl;
}



int main(int argc, char **argv)
{
  PARSE_ARGS;

  //Basic typedef's
  typedef  double PixelType;
  const unsigned int Dimension = 3;
  typedef itk::Image<PixelType, Dimension> ImageType;

  //Read in subject image
  typedef itk::ImageFileReader<ImageType> ImageReaderType;
  ImageReaderType::Pointer imageReader = ImageReaderType::New();
  imageReader->SetFileName(inputImage);
  ImageType::Pointer subject = imageReader->GetOutput();
  imageReader->Update();

  //Read in the atlas label file
  //typedef itk::Image<PixelType, Dimension> LabelAtlasType;
  //typedef itk::ImageFileReader<LabelAtlasType> LabelAtlasReaderType;
  //LabelAtlasReaderType::Pointer labelAtlasReader = LabelAtlasReaderType::New();
  //labelAtlasReader->SetFileName(labelmap);

  typedef itk::Image<unsigned char, Dimension> MaskAtlasType;

  //Read in the landmarks file
  LandmarksMapType myLandmarks = ReadSlicer3toITKLmk(landmarks);

  typedef MaskFromLandmarks<ImageType> MaskFromLandmarks;
  MaskFromLandmarks::Pointer masker = MaskFromLandmarks::New();

  masker->printHello();

  masker->SetInput(subject);
  masker->SetLandmarksFileName(landmarks);

  //WriteSmartImage<MaskAtlasType>("/scratch/aleinoff/temp/maskTest1.nii.gz", masker->GetOutput());



  //Turn Label map into binary image. Use a threshold Image filter?? or brainscut?
  //Write a new filter for this??
  //typedef itk::Image<unsigned char, Dimension> MaskAtlasType;

  /*use landmark mask
  typedef itk::BinaryThresholdImageFilter< LabelAtlasType, MaskAtlasType>  MaskFilterType;
  MaskFilterType::Pointer maskFilter = MaskFilterType::New();

  maskFilter->SetInput( labelAtlasReader->GetOutput() );
  maskFilter->SetOutsideValue(1);
  maskFilter->SetInsideValue(0);
  maskFilter->SetLowerThreshold(0);
  maskFilter->SetUpperThreshold(0);
   */


  //Write to a file
  MaskAtlasType::Pointer maskAtlas = masker->GetOutput();
  masker->Update();
  WriteImage<MaskAtlasType>(outputMask, maskAtlas);
  //Get a distance map to the Brain region:
  //typedef itk::DanielssonDistanceMapImageFilter<MaskAtlasType, ImageType, ImageType> DistanceMapFilter;
  typedef itk::SignedMaurerDistanceMapImageFilter<MaskAtlasType, ImageType> DistanceMapFilter;
  DistanceMapFilter::Pointer distanceMapFilter = DistanceMapFilter::New();
  distanceMapFilter->SetInput(maskAtlas);
  distanceMapFilter->SetSquaredDistance(false);

  //make the distance map unsigned:
  typedef itk::ThresholdImageFilter<ImageType> ThresholdFilterType;
  ThresholdFilterType::Pointer distanceThreshold = ThresholdFilterType::New();
  distanceThreshold->SetInput(distanceMapFilter->GetOutput());
  distanceThreshold->SetLower(0.0);
  distanceThreshold->SetUpper(4096);  //TODO: This should be changed to the max pixel value for the image type??? or will we always be using double for calculations??
  distanceThreshold->SetOutsideValue(0.0);



  //Write the distance map to a file so we can see what it did:
  WriteImage(distanceMapFileName, distanceThreshold->GetOutput());

  //Try to scale distance map
  typedef itk::MultiplyImageFilter<ImageType, ImageType, ImageType> ScalingFilterType;
  ScalingFilterType::Pointer distanceMapScaler = ScalingFilterType::New();
  distanceMapScaler->SetInput(distanceThreshold->GetOutput());
  distanceMapScaler->SetConstant(scaleDistanceMap);

  //Perform some kind of BSpline on Image
  const int BSplineOrder = 3;

  typedef CreateRandomBSpline<ImageType, PixelType, Dimension, BSplineOrder> BSplineCreator; //, BSTransformType> Test;
  BSplineCreator::Pointer bSplineCreator = BSplineCreator::New();
  bSplineCreator->SetInput(subject);
  bSplineCreator->SetBSplineControlPoints(bsplineControlPoints);
  bSplineCreator->SetRandMax(maxRandom);
  bSplineCreator->SetRandMin(minRandom);
  bSplineCreator->SetRandScale(scaleRandom);
  bSplineCreator->Update();

  typedef itk::BSplineTransform<PixelType, Dimension, BSplineOrder> BSTransformType;
  BSTransformType::Pointer bSpline = bSplineCreator->GetBSplineOutput();

  WriteTransform(bSplineFileName, bSpline);
//return 0;
  typedef itk::Vector<PixelType, Dimension > VectorPixelType;
  typedef itk::Image< VectorPixelType, Dimension> DisplacementFieldImageType;

  typedef CombineBSplineWithDisplacement<ImageType, DisplacementFieldImageType, PixelType, 3,3> CombinerType;

  CombinerType::Pointer combiner = CombinerType::New();
  combiner->SetBSplineInput(bSpline);
  combiner->SetInput(subject);
  combiner->SetDistanceMap(distanceMapScaler->GetOutput());
  combiner->Update();

  //write the new displacement image
  DisplacementFieldImageType* composedDisplacementField_rawPtr = combiner->GetComposedImage();
  WriteImage(smoothDisplacementName, composedDisplacementField_rawPtr);


#if 0
  //write composed displacement field into a displacement transform
  std::cout<<"printing composed image info" <<std::endl;
  combiner->GetComposedImage()->Print(std::cerr,5);
  std::cout<<"done printing composed image" <<std::endl;

  std::cout<<"printing distance map info" <<std::endl;
  combiner->GetDistanceMap()->Print(std::cerr, 0);
  std::cout<<"done printing distance map" <<std::endl;
#endif

  typedef itk::DisplacementFieldTransform<PixelType, Dimension> FinalTransformType;
  FinalTransformType::Pointer finalTransform = FinalTransformType::New();
  finalTransform->SetDisplacementField(composedDisplacementField_rawPtr);
  finalTransform->Print(std::cerr,5);

  WriteTransform(finalTransformFileName, finalTransform);

  // Apply transform to image with resampler:
  typedef itk::ResampleImageFilter<ImageType, ImageType> ResampleFilterType;
  ResampleFilterType::Pointer resampler = ResampleFilterType::New();

  typedef itk::LinearInterpolateImageFunction<ImageType, PixelType > InterpolatorType;
  InterpolatorType::Pointer interpolater = InterpolatorType::New();

  ImageType::RegionType subjectRegion = subject->GetBufferedRegion();

  resampler->SetInterpolator(interpolater);
  resampler->SetOutputSpacing(subject->GetSpacing());
  resampler->SetOutputOrigin(subject->GetOrigin());
  resampler->SetOutputDirection(subject->GetDirection());
  resampler->SetSize(subjectRegion.GetSize());
  resampler->SetOutputStartIndex(subjectRegion.GetIndex());

  resampler->SetInput(imageReader->GetOutput());
  resampler->SetTransform(finalTransform);

  WriteImage(deformedImageName, resampler->GetOutput());

  //Get the difference image
  typedef itk::SubtractImageFilter<ImageType, ImageType> SubtractFilter;
  SubtractFilter::Pointer subtractFilter = SubtractFilter::New();
  subtractFilter->SetInput1(subject);
  subtractFilter->SetInput2(resampler->GetOutput());

  //write the difference Image
  WriteImage( diffImageName, subtractFilter->GetOutput());

  std::cout << "done" << std::endl;

  return EXIT_SUCCESS;
}
