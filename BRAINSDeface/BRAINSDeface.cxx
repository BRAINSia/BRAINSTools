//
// Created by Leinoff, Alexander on 5/20/16.
//

#include <itkTransformToDisplacementFieldFilter.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <BRAINSDefaceCLP.h>
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

#include "CreateRandomBSpline.h"
#include "CombineBSplineWithDisplacement.h"

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
  typedef itk::Image<PixelType, Dimension> LabelAtlasType;
  typedef itk::ImageFileReader<LabelAtlasType> LabelAtlasReaderType;
  LabelAtlasReaderType::Pointer labelAtlasReader = LabelAtlasReaderType::New();
  labelAtlasReader->SetFileName(labelmap);

  //Read in the landmarks file
  LandmarksMapType myLandmarks = ReadSlicer3toITKLmk(landmarks);

  //Turn Label map into binary image. Use a threshold Image filter?? or brainscut?
  //Write a new filter for this??
  typedef itk::Image<unsigned char, Dimension> MaskAtlasType;
  typedef itk::BinaryThresholdImageFilter< LabelAtlasType, MaskAtlasType>  MaskFilterType;
  MaskFilterType::Pointer maskFilter = MaskFilterType::New();

  maskFilter->SetInput( labelAtlasReader->GetOutput() );
  maskFilter->SetOutsideValue(1);
  maskFilter->SetInsideValue(0);
  maskFilter->SetLowerThreshold(0);
  maskFilter->SetUpperThreshold(1);

  //Write to a file
  WriteImage(outputMask, maskFilter->GetOutput());
  //Get a distance map to the Brain region:
  typedef itk::DanielssonDistanceMapImageFilter<MaskAtlasType, ImageType, ImageType> DistanceMapFilter;
  DistanceMapFilter::Pointer distanceMapFilter = DistanceMapFilter::New();
  distanceMapFilter->SetInput(maskFilter->GetOutput());
  distanceMapFilter->InputIsBinaryOn();
  distanceMapFilter->SetSquaredDistance(false);

  //Write the distance map to a file so we can see what it did:
  WriteImage(distanceMapFileName, distanceMapFilter->GetOutput());

  //Try to scale distance map
  typedef itk::MultiplyImageFilter<ImageType, ImageType, ImageType> ScalingFilterType;
  ScalingFilterType::Pointer distanceMapScaler = ScalingFilterType::New();
  distanceMapScaler->SetInput(distanceMapFilter->GetOutput());
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


  //write composed displacement field into a displacement transform

  typedef itk::DisplacementFieldTransform<PixelType, Dimension> FinalTransformType;
  FinalTransformType::Pointer finalTransform = FinalTransformType::New();
  finalTransform->SetDisplacementField(combiner->GetComposedImage());

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
