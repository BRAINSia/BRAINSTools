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
#include <itkSignedMaurerDistanceMapImageFilter.h>
#include <itkThresholdImageFilter.h>
#include <itkTransformFileReader.h>
#include <itkCompositeTransform.h>
#include <itkTransformFactory.h>

#include "CreateRandomBSpline.h"
#include "CombineBSplineWithDisplacement.h"
#include "MaskFromLandmarksFilter.h"
#include "BRAINSRefacerUtilityFunctions.hxx"

#include "BRAINSRefacerUtilityFunctions.hxx"  //Why does this not need to be included??


int main(int argc, char **argv)
{
  PARSE_ARGS;

  if(debug_Refacer) verbose_Refacer = true;  //debug should always be verbose

  //Basic typedef's
  typedef double ProcessPixelType;
  const unsigned int Dimension = 3;
  typedef itk::Image<ProcessPixelType, Dimension> ProcessImageType;

  //Read in subject image
  typedef itk::ImageFileReader<ProcessImageType> ImageReaderType;
  ImageReaderType::Pointer imageReader = ImageReaderType::New();
  imageReader->SetFileName(inputImage);
  ProcessImageType::Pointer subject = imageReader->GetOutput();
  imageReader->Update();

  //Read the original image format

  itk::ImageIOBase::Pointer imageReaderIOBase = imageReader->GetImageIO();
  imageReaderIOBase->ReadImageInformation();
  // Note that in ImageIOBase pixel type refers to vector/scalar
  // component type refers to INT, LONG, FLOAT, etc.
  typedef itk::ImageIOBase::IOComponentType IOComponentType;
  const IOComponentType originalComponentType_ENUM = imageReaderIOBase->GetComponentType();


  typedef itk::Image<unsigned char, Dimension> ImageMaskType;
  ImageMaskType::Pointer brainMask = ImageMaskType::New();

  //Read in the atlas label file
  typedef itk::Image<ProcessPixelType, Dimension> LabelAtlasType;
  typedef itk::ImageFileReader<LabelAtlasType> LabelAtlasReaderType;
  LabelAtlasReaderType::Pointer labelAtlasReader = LabelAtlasReaderType::New();
  typedef itk::BinaryThresholdImageFilter<LabelAtlasType, ImageMaskType> MaskFilterType;
  MaskFilterType::Pointer maskFilter = MaskFilterType::New();
  LabelAtlasType::Pointer labelAtlasReaderOutput = LabelAtlasType::New();

  if (labelmapSwitch == false)
    {

    //Read in the landmarks file
    LandmarksMapType myLandmarks = ReadSlicer3toITKLmk(landmarks);

    typedef MaskFromLandmarksFilter<ProcessImageType, ImageMaskType> MaskFromLandmarksFilterType;
    MaskFromLandmarksFilterType::Pointer masker = MaskFromLandmarksFilterType::New();
    std::cout << "Generating mask from landmarks ..." <<std::endl;
    masker->SetInput(subject);
    masker->SetDebug(debug_Refacer);
    masker->SetVerbose(verbose_Refacer);
    masker->SetLandmarksFileName(landmarks);

    //Write to a file
    brainMask = masker->GetOutput();
    masker->Update();
    }
  else
    {

    std::cout << "Generating mask from labelmap" << std::endl;
    labelAtlasReader->SetFileName(labelmap);
    labelAtlasReaderOutput = labelAtlasReader->GetOutput();
    labelAtlasReader->Update();

    //resample LabelImage
    typedef itk::NearestNeighborInterpolateImageFunction<LabelAtlasType, double> NN_InterpolatorType;
    NN_InterpolatorType::Pointer NN_interpolator = NN_InterpolatorType::New();

    typedef itk::IdentityTransform<double, Dimension> IdentityTransformType;
    IdentityTransformType::Pointer identityTransform = IdentityTransformType::New();


    typedef itk::ResampleImageFilter<LabelAtlasType, LabelAtlasType> maskResamplerType;
    maskResamplerType::Pointer maskResampler = maskResamplerType::New();

    std::cout << "Resampling atlas map:" << std::endl;
    maskResampler->SetInput(labelAtlasReader->GetOutput());
    maskResampler->SetInterpolator(NN_interpolator);
    maskResampler->SetTransform(identityTransform);
    maskResampler->SetReferenceImage(subject);
    maskResampler->UseReferenceImageOn();
    maskResampler->Update();

    maskFilter->SetInput(maskResampler->GetOutput());
    maskFilter->SetOutsideValue(1);
    maskFilter->SetInsideValue(0);
    maskFilter->SetLowerThreshold(0);
    maskFilter->SetUpperThreshold(0);


    brainMask = maskFilter->GetOutput();
    brainMask->Update();
    }
  if (debug_Refacer)
    {
    WriteImage<ImageMaskType>(outputMask, brainMask);
    }
  //Get a distance map to the Brain region:
  typedef itk::SignedMaurerDistanceMapImageFilter<ImageMaskType, ProcessImageType> DistanceMapFilter;

  DistanceMapFilter::Pointer distanceMapFilter = DistanceMapFilter::New();
  std::cout << "Calculating distance map ..." << std::endl;
  distanceMapFilter->SetInput(brainMask);
  distanceMapFilter->SetSquaredDistance(false);
  ProcessImageType::Pointer myDistanceMapFilterImage = distanceMapFilter->GetOutput();
  distanceMapFilter->Update();

  //make the distance map unsigned:
  typedef itk::ThresholdImageFilter<ProcessImageType> ThresholdFilterType;
  ThresholdFilterType::Pointer distanceThreshold = ThresholdFilterType::New();
  distanceThreshold->SetInput(distanceMapFilter->GetOutput());
  distanceThreshold->SetLower(0.0);
  distanceThreshold->SetUpper(4096);  //TODO: This should be changed to the max pixel value for the image type??? or will we always be using double for calculations??
  distanceThreshold->SetOutsideValue(0.0);

  //Write the distance map to a file so we can see what it did:
  if(debug_Refacer)
    {
    WriteImage(distanceMapFileName, distanceThreshold->GetOutput());
    }
  ProcessImageType::Pointer myDistanceMapPreScaled = distanceThreshold->GetOutput();
  distanceThreshold->Update();

  //Try to scale distance map
  typedef itk::MultiplyImageFilter<ProcessImageType, ProcessImageType, ProcessImageType> ScalingFilterType;
  ScalingFilterType::Pointer distanceMapScaler = ScalingFilterType::New();
  std::cout << "Scaling distance map ..." <<std::endl;
  distanceMapScaler->SetInput(myDistanceMapPreScaled);
  distanceMapScaler->SetConstant(scaleDistanceMap);

  ProcessImageType::Pointer scaledDistanceMap = distanceMapScaler->GetOutput();
  distanceMapScaler->Update();

  //Perform some kind of BSpline on Image
  const int BSplineOrder = 3;

  typedef CreateRandomBSpline<ProcessImageType, ProcessPixelType, Dimension, BSplineOrder> BSplineCreator; //, BSTransformType> Test;
  BSplineCreator::Pointer bSplineCreator = BSplineCreator::New();
  typedef itk::BSplineTransform<ProcessPixelType, Dimension, BSplineOrder> BSTransformType;
  BSTransformType::Pointer bSpline = BSTransformType::New();

  //Stuff for reading in bspline transform
  typedef itk::TransformFileReaderTemplate< double > TransformReaderType;
  TransformReaderType::Pointer transformReader = TransformReaderType::New();


  if(!reuseBSplineSwitch)
    {
    std::cout << "Generating brand new random BSPline" << std::endl;
    bSplineCreator->SetDebug(debug_Refacer);
    bSplineCreator->SetVerbose(verbose_Refacer);
    bSplineCreator->SetInput(subject);
    bSplineCreator->SetBSplineControlPoints(bsplineControlPoints);
    bSplineCreator->SetRandMax(maxRandom);
    bSplineCreator->SetRandMin(minRandom);
    bSplineCreator->Update();
    bSpline = bSplineCreator->GetBSplineOutput();
    if(debug_Refacer || saveTransform )
      {
      WriteTransform(bSplineFileName, bSpline);
      }
    }
  else if (reuseBSplineSwitch)
    {
    std::cout << "Reusing BSpline" << std::endl;

    transformReader->SetFileName(previousBSplineFileName);

    //try catch for ioreader
    try
      {
      transformReader->Update();
      }
    catch(itk::ExceptionObject & exception)
      {
      std::cerr << "Error while reading the transform file" << std::endl;
      std::cerr << exception <<std::endl;
      return EXIT_FAILURE;
      }
    const TransformReaderType::TransformListType * transforms = transformReader->GetTransformList();
    typedef itk::CompositeTransform<double, 3> ReadCompositeTransformType;
    TransformReaderType::TransformListType::const_iterator comp_it = transforms->begin();
    if( strcmp((*comp_it)->GetNameOfClass(), "BSplineTransform") != 0 )
      {
      std::cerr << "Invalid transform given" << std::endl;
      std::cerr << "Transform type given was: " << std::endl;
      std::cerr << (*comp_it)->GetNameOfClass() << std::endl;
      std::cerr << "You should only supply BSplineTransforms" << std::endl;
      return EXIT_FAILURE;
      }

    ReadCompositeTransformType::Pointer compositeRead = static_cast<ReadCompositeTransformType *> ( (*comp_it).GetPointer() );

    // create a new bspline with the params from the composite that we read in.
    bSpline->SetFixedParameters(compositeRead->GetFixedParameters());
    bSpline->SetParameters(compositeRead->GetParameters());
    }


  typedef itk::Vector<ProcessPixelType, Dimension > VectorPixelType;
  typedef itk::Image< VectorPixelType, Dimension> DisplacementFieldProcessImageType;

  typedef CombineBSplineWithDisplacement<ProcessImageType, DisplacementFieldProcessImageType, ProcessPixelType, 3,3> CombinerType;

  CombinerType::Pointer combiner = CombinerType::New();

  std::cout << "Combining bspline with displacement ..." << std::endl;

  combiner->SetDebug(debug_Refacer);
  combiner->SetVerbose(verbose_Refacer);
  combiner->SetBSplineInput(bSpline);
  combiner->SetInput(subject);
  combiner->SetDistanceMap(distanceMapScaler->GetOutput());
  combiner->Update();

  //write the new displacement image
  DisplacementFieldProcessImageType* composedDisplacementField_rawPtr = combiner->GetComposedImage();
  if( debug_Refacer )
    {
    WriteImage(smoothDisplacementName, composedDisplacementField_rawPtr);
    }
  typedef itk::DisplacementFieldTransform<ProcessPixelType, Dimension> FinalTransformType;
  FinalTransformType::Pointer finalTransform = FinalTransformType::New();
  finalTransform->SetDisplacementField(composedDisplacementField_rawPtr);

  if( debug_Refacer )
    {
    WriteTransform(finalTransformFileName, finalTransform);
    }

  // Apply transform to image with resampler:
  typedef itk::ResampleImageFilter<ProcessImageType, ProcessImageType> ResampleFilterType;
  ResampleFilterType::Pointer resampler = ResampleFilterType::New();

  typedef itk::LinearInterpolateImageFunction<ProcessImageType, ProcessPixelType > InterpolatorType;
  InterpolatorType::Pointer interpolater = InterpolatorType::New();

  ProcessImageType::RegionType subjectRegion = subject->GetBufferedRegion();

  std::cout << "Refacing image ..." << std::endl;

  resampler->SetInterpolator(interpolater);
  resampler->SetReferenceImage(imageReader->GetOutput());
  resampler->UseReferenceImageOn();
  resampler->SetOutputStartIndex(subjectRegion.GetIndex());

  resampler->SetInput(imageReader->GetOutput());
  resampler->SetTransform(finalTransform);

  //WriteImage(deformedImageName, resampler->GetOutput());
  ConvertAndSave<ProcessImageType, Dimension>( deformedImageName, resampler->GetOutput(), originalComponentType_ENUM);

  //write the difference Image
  if( debug_Refacer )
    {
    //Get the difference image
    typedef itk::SubtractImageFilter<ProcessImageType, ProcessImageType> SubtractFilter;
    SubtractFilter::Pointer subtractFilter = SubtractFilter::New();
    subtractFilter->SetInput1(subject);
    subtractFilter->SetInput2(resampler->GetOutput());
    WriteImage( diffImageName, subtractFilter->GetOutput());
    }

  std::cout << "done" << std::endl;

  return EXIT_SUCCESS;
}
