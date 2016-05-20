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

inline double myRandom()
{
  const int min = -5;
  const int max = 5;
  const int range = max - min;
  return static_cast< double >( rand() % range +1 - max )*0.05;
}

template< typename TImageType >
void writeAnImage(std::string filename, TImageType* image)
{

  typedef itk::ImageFileWriter<TImageType> FileWriterType;
  typename FileWriterType::Pointer fileWriter = FileWriterType::New();
  fileWriter->SetInput(image);
  fileWriter->SetFileName(filename);
  fileWriter->Update();
}

int main(int argc, char **argv)
{
  PARSE_ARGS;

  //Read in the imagefile. Right now assumed to be t1 weighted image

  typedef  double PixelType;
  //TODO: how to activate cpp11 in BRAINSTools so I can use constexpr???
  const unsigned int Dimension = 3;

  typedef itk::Image<PixelType, Dimension> ImageType;
  typedef itk::ImageFileReader<ImageType> ImageReaderType;
  ImageReaderType::Pointer imageReader = ImageReaderType::New();
  imageReader->SetFileName(inputImage);

  //Extract subject image from reader
  ImageType::Pointer subject = imageReader->GetOutput();
  imageReader->Update();
  //


  //Read in the atlas label file
  typedef itk::Image<PixelType, Dimension> LabelAtlasType;
  typedef itk::ImageFileReader<LabelAtlasType> LabelAtlasReaderType;
  LabelAtlasReaderType::Pointer labelAtlasReader = LabelAtlasReaderType::New();
  labelAtlasReader->SetFileName(labelmap);

  //Read in the landmarks file
  LandmarksMapType myLandmarks = ReadSlicer3toITKLmk(landmarks);


  //Turn Label map into binary image. Use a threshold Image filter?? or brainscut?
  typedef itk::Image<unsigned char, Dimension> MaskAtlasType;
  typedef itk::BinaryThresholdImageFilter< LabelAtlasType, MaskAtlasType>  MaskFilterType;
  MaskFilterType::Pointer maskFilter = MaskFilterType::New();

  maskFilter->SetInput( labelAtlasReader->GetOutput() );
  maskFilter->SetOutsideValue(1);
  maskFilter->SetInsideValue(0);
  maskFilter->SetLowerThreshold(0);
  maskFilter->SetUpperThreshold(1);

  //Write to a file
  writeAnImage(outputMask, maskFilter->GetOutput());

  //Get a distance map to the Brain region:
  //TODO: This should be changed to the other kind of distance map that is faster(Mauer I think it is called???)
  typedef itk::DanielssonDistanceMapImageFilter<MaskAtlasType, ImageType, ImageType> DistanceMapFilter;
  DistanceMapFilter::Pointer distanceMapFilter = DistanceMapFilter::New();
  distanceMapFilter->SetInput(maskFilter->GetOutput());
  distanceMapFilter->InputIsBinaryOn();


  //Write the distance map to a file so we can see what it did:
  writeAnImage(distanceMapFileName, distanceMapFilter->GetOutput());


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

  //Setup BSpline
  const unsigned int BSplineOrder = 3;
  typedef itk::BSplineTransform<PixelType, Dimension, BSplineOrder> BSplineTransform;

  BSplineTransform::Pointer bSpline = BSplineTransform::New();
#if 0 //Set this up later
  typedef itk::BSplineTransformInitializer<BSplineTransform> BSplineInitType;


  BSplineInitType::Pointer bsi = BSplineInitType::New();
  bsi->SetImage(subject);
  bsi->SetTransform(bSpline);
  bsi->Update();
#endif
  //Set BSpline basic parameters
  const unsigned int BSplineControlPoints = 8;
  BSplineTransform::MeshSizeType meshSize;                  //Setup a mesh that contains the number of controlpoints
  meshSize.Fill(BSplineControlPoints-BSplineOrder);         //Inspired from itk example "BSplineWarping2.cxx"

  bSpline->SetTransformDomainMeshSize(meshSize);            //TODO: ask if it is possible to have "uneven" control points.
  // EG. 8 on LR axis 7 on SI 6 on AP. We did this in simple itk
  //so it should be possible in cpp itk


#if 1 // we should change this so we sue BSplineTransform initializer instead
  typedef ImageType::RegionType ImageRegionType;
//  ImageRegionType subjectRegion = subject->GetLargestPossibleRegion();
  ImageRegionType subjectRegion = subject->GetBufferedRegion();

  bSpline->SetTransformDomainOrigin(subject->GetOrigin());           //Origin
  bSpline->SetTransformDomainDirection(subject->GetDirection());     //Direction
  bSpline->SetTransformDomainPhysicalDimensions((                    //PhysicalDimensions
    subject->GetSpacing()[0]*(subjectRegion.GetSize()[0]-1),         //Should all be set to the same as the subject Image
    subject->GetSpacing()[1]*(subjectRegion.GetSize()[1]-1),
    subject->GetSpacing()[2]*(subjectRegion.GetSize()[2]-1)
    ));
#endif

  //Get the number of paramaters/nodes required for this BSpline
  const unsigned int numberOfParameters = bSpline->GetNumberOfParameters();
  const unsigned int numberOfNodes = numberOfParameters/Dimension;

  //print out
  std::cout << "Number of params: " << numberOfParameters << std::endl;
  std::cout << "Number of nodes:  " << numberOfNodes << std::endl;

  //Setup a paramaters variable for the bspline
  BSplineTransform::ParametersType bSplineParams( numberOfParameters );

  //  From ITK Example "BSplineWarping2"
  //  The B-spline grid should now be fed with coeficients at each node. Since
  //  this is a two dimensional grid, each node should receive two coefficients.
  //  Each coefficient pair is representing a displacement vector at this node.
  //  The coefficients can be passed to the B-spline in the form of an array where
  //  the first set of elements are the first component of the displacements for
  //  all the nodes, and the second set of elemets is formed by the second
  //  component of the displacements for all the nodes.

  //  In the ITK Example, the read the points in from a file. Here, they will be
  //  generated randomly. This should put the xyz coordinates in the correct space
  //  a better way would probably be to use the image coefficient array of the bspline,
  //  but for now I will use the method from the ITK example

  std::srand(time(nullptr));

  for( unsigned int n = 0; n < numberOfNodes; ++ n)
    {
    bSplineParams[n] = myRandom();
    bSplineParams[n + numberOfNodes] = myRandom();      // "y" coord;
    bSplineParams[n + numberOfNodes * 2] = myRandom();  // "z" coord;
                                                                  // TODO: x,y,z seem like they are the wrong coordinate system. Get a better model
    }

  bSpline->SetParameters(bSplineParams);

//Write the bspline transform
  typedef itk::TransformFileWriter TransformWriterType;
  TransformWriterType::Pointer transformWriter = TransformWriterType::New();
  transformWriter->SetInput(bSpline);
  transformWriter->SetFileName(bSplineFileName);
  transformWriter->Update();

  std::cout << "Printing bSpline paramaters" << std::endl;
  std::cout << bSpline->GetParameters() << std::endl;
  std::cout <<"Printing bSpline info" << std::endl;
  bSpline->Print(std::cout,0);
  std::cout << std::endl <<std::endl<<std::endl;


  //Get the displacement field from the bspline transform



  //The displacement field is in a vector image
  typedef itk::Vector<PixelType, Dimension > VectorPixelType;
  typedef itk::Image< VectorPixelType, Dimension> DisplacementFieldImageType;

  typedef itk::TransformToDisplacementFieldFilter< DisplacementFieldImageType, PixelType> TransformToDisplacementFilterType;
  TransformToDisplacementFilterType::Pointer bSplineDisplacementFieldGenerator = TransformToDisplacementFilterType::New();
  bSplineDisplacementFieldGenerator->UseReferenceImageOn();
  bSplineDisplacementFieldGenerator->SetReferenceImage(subject);
  bSplineDisplacementFieldGenerator->SetTransform(bSpline);
  //bSplineDisplacementFieldGenerator->Print(std::cout,0);

  writeAnImage(bSplineDisplacementFileName, bSplineDisplacementFieldGenerator->GetOutput());

  std::cout<<"done writing initial bSpline displacmentField"<<std::endl;

  std::cout<<"Extracting component images from displacement field"<<std::endl;
  //multiply the displacement field by the distance map to get the "smooth displacement that doesn't affect the brain
  //first extrace scalar elemnts from vector image
  typedef itk::VectorIndexSelectionCastImageFilter<DisplacementFieldImageType, ImageType> ImageExtractionFilterType;
  ImageExtractionFilterType::Pointer xTractDisplacementFilter = ImageExtractionFilterType::New();
  ImageExtractionFilterType::Pointer yTractDisplacementFilter = ImageExtractionFilterType::New();
  ImageExtractionFilterType::Pointer zTractDisplacementFilter = ImageExtractionFilterType::New();

  xTractDisplacementFilter->SetInput(bSplineDisplacementFieldGenerator->GetOutput());
  xTractDisplacementFilter->SetIndex(0);
  ImageType::Pointer   xDisplacement = xTractDisplacementFilter->GetOutput();


  //writeAnImage("/scratch/aleinoff/defaceOutput/testImageVectorIndexSelectionCastFilterOutput.nii.gz", xTractDisplacementFilter->GetOutput());
  //std::cout<<"done writing test image"<<std::endl;
  //return 0;


  yTractDisplacementFilter->SetInput(bSplineDisplacementFieldGenerator->GetOutput());
  yTractDisplacementFilter->SetIndex(1);
  ImageType::Pointer  yDisplacement = yTractDisplacementFilter->GetOutput();

  zTractDisplacementFilter->SetInput(bSplineDisplacementFieldGenerator->GetOutput());
  zTractDisplacementFilter->SetIndex(2);
  ImageType::Pointer  zDisplacement = zTractDisplacementFilter->GetOutput();

  //multiply by distancemap

  std::cout<<"Multiplying bSplineDisplacement by Distance map"<<std::endl;

  typedef itk::MultiplyImageFilter<ImageType, ImageType, ImageType> MultiplyFilterType;
  MultiplyFilterType::Pointer xMult = MultiplyFilterType::New();
  MultiplyFilterType::Pointer yMult = MultiplyFilterType::New();
  MultiplyFilterType::Pointer zMult = MultiplyFilterType::New();

  xMult->SetInput1(xDisplacement);
  xMult->SetInput2(distanceMapFilter->GetOutput());
  ImageType::Pointer xMultImage = xMult->GetOutput();
  xMult->Update();


  yMult->SetInput1(yDisplacement);
  yMult->SetInput2(distanceMapFilter->GetOutput());
  ImageType::Pointer yMultImage = yMult->GetOutput();
  yMult->Update();


  zMult->SetInput1(zDisplacement);
  zMult->SetInput2(distanceMapFilter->GetOutput());
  ImageType::Pointer zMultImage = zMult->GetOutput();
  zMult->Update();

  std::cout<<"Composing new image from displacement and bSpline product components"<<std::endl;

  typedef itk::ComposeImageFilter<ImageType, DisplacementFieldImageType> ComposeFilterType;
  ComposeFilterType::Pointer composeDisplacements = ComposeFilterType::New();
  composeDisplacements->SetInput(0, xMult->GetOutput());
  composeDisplacements->SetInput(1, yMult->GetOutput());
  composeDisplacements->SetInput(2, zMult->GetOutput());

  //write the new displacement image
  writeAnImage(smoothDisplacementName, composeDisplacements->GetOutput());

  DisplacementFieldImageType::Pointer composedDisplacementField = composeDisplacements->GetOutput();
  composedDisplacementField->Update();

  //write composed displacement field into a displacement transform

  typedef itk::DisplacementFieldTransform<PixelType, Dimension> FinalTransformType;
  FinalTransformType::Pointer finalTransform = FinalTransformType::New();
  finalTransform->SetDisplacementField(composedDisplacementField);




  //DisplacementFieldImageType::Pointer bSplineDisplacementField = bSplineDisplacementFieldGenerator->GetOutput();

  //write the displacement field


  // Apply transform to image with resampler:
  typedef itk::ResampleImageFilter<ImageType, ImageType> ResampleFilterType;
  ResampleFilterType::Pointer resampler = ResampleFilterType::New();

  typedef itk::LinearInterpolateImageFunction<ImageType, PixelType > InterpolatorType;
  InterpolatorType::Pointer interpolater = InterpolatorType::New();

  resampler->SetInterpolator(interpolater);
  resampler->SetOutputSpacing(subject->GetSpacing());
  resampler->SetOutputOrigin(subject->GetOrigin());
  resampler->SetOutputDirection(subject->GetDirection());
  resampler->SetSize(subjectRegion.GetSize());
  resampler->SetOutputStartIndex(subjectRegion.GetIndex());

  resampler->SetInput(imageReader->GetOutput());
  //resampler->SetTransform(bSpline);
  resampler->SetTransform(finalTransform);



  writeAnImage(deformedImageName,resampler->GetOutput());

  //Get the difference image
  typedef itk::SubtractImageFilter<ImageType, ImageType> SubtractFilter;
  SubtractFilter::Pointer subtractFilter = SubtractFilter::New();
  subtractFilter->SetInput1(subject);
  subtractFilter->SetInput2(resampler->GetOutput());

  //write the difference Image
  writeAnImage("/scratch/aleinoff/defaceOutput/diffImage.nii.gz", subtractFilter->GetOutput());

  /*
  typedef itk::ImageFileWriter<ImageType> DeformedSubjectFileWriterType;
  DeformedSubjectFileWriterType::Pointer deformedWriter = DeformedSubjectFileWriterType::New();
  deformedWriter->SetFileName(deformedImageName);
  deformedWriter->SetInput(resampler->GetOutput());
  deformedWriter->Update();
   */
  std::cout << "Finished writing file " << std::endl;
  std::cout << "done" << std::endl;


  return EXIT_SUCCESS;
}