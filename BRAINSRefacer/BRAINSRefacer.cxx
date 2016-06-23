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


#include "itkCastImageFilter.h"
#include "itkRoundImageFilter.h"

//Convenience function to write images
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

//Convenience function to write images
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
// One ugly function
template< typename TInputImageType, int TDimension>
int ConvertAndSave(std::string fileName, typename TInputImageType::Pointer image, itk::ImageIOBase::IOComponentType InputComponentType);



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
  typedef  double                                                                 ProcessPixelType;
  const unsigned int                                                              Dimension = 3;
  typedef itk::Image<ProcessPixelType, Dimension>                                 ProcessImageType;

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
  typedef itk::ImageIOBase::IOComponentType                                       IOComponentType;
  const IOComponentType originalComponentType_ENUM = imageReaderIOBase->GetComponentType();

  //Read in the atlas label file
  //typedef itk::Image<ProcessPixelType, Dimension> LabelAtlasType;
  //typedef itk::ImageFileReader<LabelAtlasType> LabelAtlasReaderType;
  //LabelAtlasReaderType::Pointer labelAtlasReader = LabelAtlasReaderType::New();
  //labelAtlasReader->SetFileName(labelmap);

  typedef itk::Image<unsigned char, Dimension> MaskAtlasType;

  //Read in the landmarks file
  LandmarksMapType myLandmarks = ReadSlicer3toITKLmk(landmarks);

  typedef MaskFromLandmarks<ProcessImageType> MaskFromLandmarks;
  MaskFromLandmarks::Pointer masker = MaskFromLandmarks::New();

  masker->printHello();

  masker->SetInput(subject);
  masker->SetLandmarksFileName(landmarks);

  //Write to a file
  MaskAtlasType::Pointer maskAtlas = masker->GetOutput();
  masker->Update();
  WriteImage<MaskAtlasType>(outputMask, maskAtlas);
  //Get a distance map to the Brain region:
  typedef itk::SignedMaurerDistanceMapImageFilter<MaskAtlasType, ProcessImageType> DistanceMapFilter;

  DistanceMapFilter::Pointer distanceMapFilter = DistanceMapFilter::New();
  distanceMapFilter->SetInput(maskAtlas);
  distanceMapFilter->SetSquaredDistance(false);

  //make the distance map unsigned:
  typedef itk::ThresholdImageFilter<ProcessImageType> ThresholdFilterType;
  ThresholdFilterType::Pointer distanceThreshold = ThresholdFilterType::New();
  distanceThreshold->SetInput(distanceMapFilter->GetOutput());
  distanceThreshold->SetLower(0.0);
  distanceThreshold->SetUpper(4096);  //TODO: This should be changed to the max pixel value for the image type??? or will we always be using double for calculations??
  distanceThreshold->SetOutsideValue(0.0);



  //Write the distance map to a file so we can see what it did:
  WriteImage(distanceMapFileName, distanceThreshold->GetOutput());

  //Try to scale distance map
  typedef itk::MultiplyImageFilter<ProcessImageType, ProcessImageType, ProcessImageType> ScalingFilterType;
  ScalingFilterType::Pointer distanceMapScaler = ScalingFilterType::New();
  distanceMapScaler->SetInput(distanceThreshold->GetOutput());
  distanceMapScaler->SetConstant(scaleDistanceMap);

  //Perform some kind of BSpline on Image
  const int BSplineOrder = 3;

  typedef CreateRandomBSpline<ProcessImageType, ProcessPixelType, Dimension, BSplineOrder> BSplineCreator; //, BSTransformType> Test;
  BSplineCreator::Pointer bSplineCreator = BSplineCreator::New();
  bSplineCreator->SetInput(subject);
  bSplineCreator->SetBSplineControlPoints(bsplineControlPoints);
  bSplineCreator->SetRandMax(maxRandom);
  bSplineCreator->SetRandMin(minRandom);
  bSplineCreator->SetRandScale(scaleRandom);
  bSplineCreator->Update();

  typedef itk::BSplineTransform<ProcessPixelType, Dimension, BSplineOrder> BSTransformType;
  BSTransformType::Pointer bSpline = bSplineCreator->GetBSplineOutput();

  WriteTransform(bSplineFileName, bSpline);

  typedef itk::Vector<ProcessPixelType, Dimension > VectorPixelType;
  typedef itk::Image< VectorPixelType, Dimension> DisplacementFieldProcessImageType;

  typedef CombineBSplineWithDisplacement<ProcessImageType, DisplacementFieldProcessImageType, ProcessPixelType, 3,3> CombinerType;

  CombinerType::Pointer combiner = CombinerType::New();
  combiner->SetBSplineInput(bSpline);
  combiner->SetInput(subject);
  combiner->SetDistanceMap(distanceMapScaler->GetOutput());
  combiner->Update();

  //write the new displacement image
  DisplacementFieldProcessImageType* composedDisplacementField_rawPtr = combiner->GetComposedImage();
  WriteImage(smoothDisplacementName, composedDisplacementField_rawPtr);

  typedef itk::DisplacementFieldTransform<ProcessPixelType, Dimension> FinalTransformType;
  FinalTransformType::Pointer finalTransform = FinalTransformType::New();
  finalTransform->SetDisplacementField(composedDisplacementField_rawPtr);
  finalTransform->Print(std::cerr,5);

  WriteTransform(finalTransformFileName, finalTransform);

  // Apply transform to image with resampler:
  typedef itk::ResampleImageFilter<ProcessImageType, ProcessImageType> ResampleFilterType;
  ResampleFilterType::Pointer resampler = ResampleFilterType::New();

  typedef itk::LinearInterpolateImageFunction<ProcessImageType, ProcessPixelType > InterpolatorType;
  InterpolatorType::Pointer interpolater = InterpolatorType::New();

  ProcessImageType::RegionType subjectRegion = subject->GetBufferedRegion();

  resampler->SetInterpolator(interpolater);
  resampler->SetOutputSpacing(subject->GetSpacing());
  resampler->SetOutputOrigin(subject->GetOrigin());
  resampler->SetOutputDirection(subject->GetDirection());
  resampler->SetSize(subjectRegion.GetSize());
  resampler->SetOutputStartIndex(subjectRegion.GetIndex());

  resampler->SetInput(imageReader->GetOutput());
  resampler->SetTransform(finalTransform);

  //WriteImage(deformedImageName, resampler->GetOutput());
  ConvertAndSave<ProcessImageType, Dimension>( deformedImageName, resampler->GetOutput(), originalComponentType_ENUM);

  //Get the difference image
  typedef itk::SubtractImageFilter<ProcessImageType, ProcessImageType> SubtractFilter;
  SubtractFilter::Pointer subtractFilter = SubtractFilter::New();
  subtractFilter->SetInput1(subject);
  subtractFilter->SetInput2(resampler->GetOutput());

  //write the difference Image
  WriteImage( diffImageName, subtractFilter->GetOutput());

  std::cout << "done" << std::endl;

  return EXIT_SUCCESS;
}

template<typename TInputImageType, typename TOutputImageType>
void RoundAndWriteImage( std::string fileName, typename TInputImageType::Pointer image)
{
  auto rounder = itk::RoundImageFilter<TInputImageType, TOutputImageType>::New();
  rounder->SetInput(image);
  WriteSmartImage<TOutputImageType>(fileName, rounder->GetOutput());
}

template<typename TInputImageType, typename TOutputImageType>
void CastAndWriteImage( std::string fileName, typename TInputImageType::Pointer image)
{
  auto caster = itk::CastImageFilter<TInputImageType, TOutputImageType>::New();
  caster->SetInput(image);
  WriteSmartImage<TOutputImageType>(fileName, caster->GetOutput());
}

template< typename TInputImageType, int TDimension>
int ConvertAndSave(std::string fileName, typename TInputImageType::Pointer image, itk::ImageIOBase::IOComponentType OutputComponentType_ENUM)
{
  switch( OutputComponentType_ENUM )
    {

    //INTEGER like types  --- round and write
    case itk::ImageIOBase::SHORT:
      {
      typedef itk::Image<short, TDimension> OutputType;
      RoundAndWriteImage<TInputImageType, OutputType>(fileName, image);
      break;
      }

    case itk::ImageIOBase::UCHAR:
      {
      typedef itk::Image<unsigned char, TDimension> OutputType;
      RoundAndWriteImage<TInputImageType, OutputType>(fileName, image);
      break;
      }

    case itk::ImageIOBase::CHAR:
      {
      typedef itk::Image<char, TDimension> OutputType;
      RoundAndWriteImage<TInputImageType, OutputType>(fileName, image);
      break;
      }

    case itk::ImageIOBase::USHORT:
      {
      typedef itk::Image<unsigned short, TDimension> OutputType;
      RoundAndWriteImage<TInputImageType, OutputType>(fileName, image);
      break;
      }

    case itk::ImageIOBase::UINT:
      {
      typedef itk::Image<unsigned int, TDimension> OutputType;
      RoundAndWriteImage<TInputImageType, OutputType> (fileName, image);
      break;
      }

    case itk::ImageIOBase::INT:
      {
      typedef itk::Image<int, TDimension> OutputType;
      RoundAndWriteImage<TInputImageType, OutputType> (fileName, image);
      break;
      }

    case itk::ImageIOBase::ULONG:
      {
      typedef itk::Image<unsigned long, TDimension> OutputType;
      RoundAndWriteImage<TInputImageType, OutputType> (fileName, image);
      break;
      }

    case itk::ImageIOBase::LONG:
      {
      typedef itk::Image<long, TDimension> OutputType;
      RoundAndWriteImage<TInputImageType, OutputType> (fileName, image);
      break;
      }


    //FLOAT LIKE TYPES --- cast and write
    case itk::ImageIOBase::FLOAT:
      {
      typedef itk::Image<float, TDimension> OutputType;
      CastAndWriteImage<TInputImageType, OutputType> (fileName, image);
      break;
      }

    case itk::ImageIOBase::DOUBLE:
      {
      typedef itk::Image<double, TDimension> OutputType;
      CastAndWriteImage<TInputImageType, OutputType> (fileName, image);
      break;
      }


    //ERROR producing cases
    case itk::ImageIOBase::UNKNOWNCOMPONENTTYPE:
      {
      std::cerr << "Unknown component type: ..." << std::endl;
      return EXIT_FAILURE;
      }

    default:
      {
      std::cerr << "Unprocessed Case: " << itk::ImageIOBase::GetComponentTypeAsString(OutputComponentType_ENUM) << std::endl;
      return EXIT_FAILURE;
      }
    }

  return EXIT_SUCCESS;
}