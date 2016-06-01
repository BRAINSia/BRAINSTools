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
#include <stdlib.h>
#include <time.h>
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
  typedef itk::ImageFileWriter<TImageType> FileWriterType;
  typename FileWriterType::Pointer fileWriter = FileWriterType::New();

  fileWriter->SetInput(image);
  fileWriter->SetFileName(filename);
  fileWriter->Update();
}

//Convienience function to write transforms
template< typename TTransformType >
void WriteTransform(std::string transformFileName, TTransformType transform )
{
  typedef itk::TransformFileWriter TransformWriterType;
  TransformWriterType::Pointer transformWriter = TransformWriterType::New();

  transformWriter->SetInput(transform);
  transformWriter->SetFileName(transformFileName);
  transformWriter->Update();
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

  //Write the distance map to a file so we can see what it did:
  WriteImage(distanceMapFileName, distanceMapFilter->GetOutput());


/*
  //scale the distance map to reduce displacement:
  typedef itk::DivideImageFilter<ImageType, ImageType, ImageType> ScaledDistanceMapFilter;
  ScaledDistanceMapFilter::Pointer scaledDistanceMapFilter = ScaledDistanceMapFilter::New();
  scaledDistanceMapFilter->SetConstant1(10);
  scaledDistanceMapFilter->SetInput(distanceMapFilter->GetOutput());
  //write it to a file:
  DistanceMapWriterType::Pointer scaledDistanceMapWriter = DistanceMapWriterType::New();
  scaledDistanceMapWriter->SetInput(scaledDistanceMapFilter->GetOutput());
  scaledDistanceMapWriter->SetFileName(distanceMapFileName);
  scaledDistanceMapWriter->Update();


  //Get the new scaled distance map as an Image:
  ImageType::Pointer distanceMap = scaledDistanceMapFilter->GetOutput();
  scaledDistanceMapFilter->Update();
*/


  //Perform some kind of BSpline on Image
  const int BSplineOrder = 3;
  const int BSplineControlPoints = 8;

  //TODO: call function

  /*
   * template<typename TImageType, typename TBSplineType>
TBSplineType createRandomBSpline(TImageType subject, const int Dimension, const int BSplineOrder, const int BSplineControlPoints)
{
   */
  typedef itk::BSplineTransform<PixelType, Dimension, BSplineOrder> BSTransformType;

  //BSTransformType::Pointer bSpline = createRandomBSpline2<ImageType, BSTransformType>(subject, Dimension, BSplineOrder, BSplineControlPoints );
  typedef CreateRandomBSpline<ImageType, PixelType, 3, 3> Test; //, BSTransformType> Test;
  Test::Pointer myTest = Test::New();
  myTest->SetInput(subject);
  myTest->SetBSplineControlPoints(8);

  myTest->Update();
  BSTransformType::Pointer bSpline = myTest->GetBSplineOutput();

  myTest->Print(std::cerr,5);


  WriteTransform(bSplineFileName, bSpline);
  std::cout << "Printing bSpline paramaters" << std::endl;
  std::cout << bSpline->GetParameters() << std::endl;


  std::cout <<"Printing bSpline info" << std::endl;


  bSpline->Print(std::cout,0);
  std::cout<<   "printed bspline info"<<std::endl;
  std::cout << std::endl <<std::endl<<std::endl;

  //return 0;
  //Get the displacement field from the bspline transform

  typedef itk::Vector<PixelType, Dimension > VectorPixelType;
  typedef itk::Image< VectorPixelType, Dimension> DisplacementFieldImageType;

 // DisplacementFieldImageType::Pointer bSplineDistanceCombo = BSplineDisplacementMultiplier<ImageType, 3>(subject, bSpline, distanceMapFilter->GetOutput());
 // DisplacementFieldImageType * bSplineDistanceCombo_Rawptr = bSplineDistanceCombo;

  typedef CombineBSplineWithDisplacement<ImageType, DisplacementFieldImageType, PixelType, 3,3> CombinerType;
  CombinerType::Pointer combiner = CombinerType::New();
  combiner->SetBSplineInput(bSpline);
  combiner->SetInput(subject);
  combiner->SetDistanceMap(distanceMapFilter->GetOutput());
  combiner->Update();

  DisplacementFieldImageType* dispCombo = combiner->GetComposedImage();
 // WriteImage("/scratch/aleinoff/defaceOutput/test.nii.gz", dispCombo);



//refactoring took place here

  //write the new displacement image

  DisplacementFieldImageType* composedDisplacementField_rawPtr = combiner->GetComposedImage();
  WriteImage(smoothDisplacementName, composedDisplacementField_rawPtr);
  //composedDisplacementField->Update();




 // DisplacementFieldImageType::Pointer composedDisplacementField = composeDisplacements->GetOutput();
 // composedDisplacementField->Update();

  //write composed displacement field into a displacement transform

  typedef itk::DisplacementFieldTransform<PixelType, Dimension> FinalTransformType;
  FinalTransformType::Pointer finalTransform = FinalTransformType::New();
  finalTransform->SetDisplacementField(combiner->GetComposedImage());




  //DisplacementFieldImageType::Pointer bSplineDisplacementField = bSplineDisplacementFieldGenerator->GetOutput();

  //write the displacement field


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
  //resampler->SetTransform(bSpline);
  resampler->SetTransform(finalTransform);


  WriteImage(deformedImageName, resampler->GetOutput());

  //Get the difference image
  typedef itk::SubtractImageFilter<ImageType, ImageType> SubtractFilter;
  SubtractFilter::Pointer subtractFilter = SubtractFilter::New();
  subtractFilter->SetInput1(subject);
  subtractFilter->SetInput2(resampler->GetOutput());

  //write the difference Image
  WriteImage( diffImageName , subtractFilter->GetOutput());

  std::cout << "Finished writing file " << std::endl;
  std::cout << "done" << std::endl;


  return EXIT_SUCCESS;
}
