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
//
// This is a simple program for mattes mutual information metric
// evaluation over two input images and a BSpline transform.
//
//
#include "itkMattesMutualInformationImageToImageMetricv4.h"
#include "itkBSplineTransform.h"
#include "itkImageRegionConstIteratorWithIndex.h"

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkTransformFileReader.h"

int main(int argc, char* argv[])
{
  if( argc < 3 )
    {
    std::cerr << "Missing Parameters " << std::endl;
    std::cerr << "Usage: " << argv[0];
    std::cerr << " fixedImageFile  movingImageFile ";
    std::cerr << "[transformFile]" << std::endl;
    std::cerr << "[NumberOfHistogramBins]" << std::endl;
    std::cerr << "[NumberOfSamples]" << std::endl;
    return EXIT_FAILURE;
    }

  constexpr unsigned int Dimension = 3;
  typedef  float           PixelType;
  typedef itk::Image< PixelType, Dimension > FixedImageType;
  typedef itk::Image< PixelType, Dimension > MovingImageType;

  typedef itk::ImageFileReader< FixedImageType  >  FixedImageReaderType;
  FixedImageReaderType::Pointer fixedReader = FixedImageReaderType::New();
  fixedReader->SetFileName( argv[1] );
  FixedImageType::Pointer fixedImage = fixedReader->GetOutput();
  fixedReader->Update();

  typedef itk::ImageFileReader< MovingImageType  >  MovingImageReaderType;
  MovingImageReaderType::Pointer movingReader = MovingImageReaderType::New();
  movingReader->SetFileName( argv[2] );
  MovingImageType::Pointer movingImage = movingReader->GetOutput();
  movingReader->Update();

  typedef itk::BSplineTransform< double, Dimension > TransformType;
  TransformType::Pointer transform = TransformType::New();

  if( argc > 3 )
    {
    std::cout << "Read transform file from the disk ..." << std::endl;
    itk::TransformFileReader::Pointer transReader = itk::TransformFileReader::New();
    transReader->SetFileName( argv[3] );
    try
      {
      transReader->Update();
      }
    catch( itk::ExceptionObject & e )
      {
      std::cerr << "Exception caught!" << std::endl;
      std::cerr << e << std::endl;
      return EXIT_FAILURE;
      }
    transform = dynamic_cast<TransformType *>( transReader->GetTransformList()->begin()->GetPointer() );
    }
  else
    {
    std::cout << "No input transform file ... transform file is set to identity" << std::endl;
    transform->SetIdentity();
    }

  typedef itk::MattesMutualInformationImageToImageMetricv4<FixedImageType, MovingImageType>     MIMetricType;
  MIMetricType::Pointer metric = MIMetricType::New();

  metric->SetVirtualDomainFromImage( fixedImage );
  metric->SetFixedImage( fixedImage );
  metric->SetMovingImage( movingImage );
  metric->SetMovingTransform( transform );

  unsigned int NoHistBins = 50;
  if ( argc > 4)
    {
    NoHistBins = atoi( argv[4] );
    }

  metric->SetNumberOfHistogramBins( NoHistBins );

  metric->SetUseFixedImageGradientFilter( false );
  metric->SetUseMovingImageGradientFilter( false );

  //  REGULAR sampling
  std::cout << "Use regular sampling ..." << std::endl;
  unsigned long NoOfSamples = 14000;
  if ( argc > 5)
    {
    NoOfSamples = atoi( argv[5] );
    }

  const unsigned long numberOfAllSamples = fixedImage->GetBufferedRegion().GetNumberOfPixels();
  const double samplingPercentage = static_cast<double>(NoOfSamples)/numberOfAllSamples;
  std::cout << "numberOfAllSamples: " << numberOfAllSamples << std::endl;
  std::cout << "Sampling percentage: " << samplingPercentage << std::endl;

  typedef MIMetricType::FixedSampledPointSetType          MetricSamplePointSetType;
  MetricSamplePointSetType::Pointer samplePointSet = MetricSamplePointSetType::New();
  samplePointSet->Initialize();
  typedef MetricSamplePointSetType::PointType SamplePointType;
  unsigned long index = 0;

  const unsigned long sampleCount = static_cast<unsigned long>( std::ceil( 1.0 / samplingPercentage ) );
  std::cout << "sample count: " << sampleCount << std::endl;
  unsigned long count = sampleCount;
  itk::ImageRegionConstIteratorWithIndex<FixedImageType> It( fixedImage, fixedImage->GetBufferedRegion() );
  for( It.GoToBegin(); !It.IsAtEnd(); ++It )
    {
    if( count == sampleCount )
      {
      count=0; //Reset counter
      SamplePointType point;
      fixedImage->TransformIndexToPhysicalPoint( It.GetIndex(), point );
      samplePointSet->SetPoint( index, point );
      ++index;
      }
    ++count;
    }

  metric->SetFixedSampledPointSet( samplePointSet );
  metric->SetUseFixedSampledPointSet( true );
  ////////////////////////

  metric->Initialize();

  MIMetricType::MeasureType measure;
  MIMetricType::DerivativeType derivative( metric->GetTransform()->GetNumberOfParameters() );
  try
    {
    metric->GetValueAndDerivative( measure, derivative );
    }
  catch ( itk::ExceptionObject & err )
    {
    std::cerr << "Exception caught!" << std::endl;
    std::cerr << err << std::endl;
    return EXIT_FAILURE;
    }
  std::cout << "NumberOfValidPoints: " << metric->GetNumberOfValidPoints() << " of "
            << metric->GetVirtualRegion().GetNumberOfPixels() << std::endl;
  std::cout << "measure: " << measure << std::endl;
//  std::cout << "Derivatives: " << derivative << std::endl;
//  std::cout << "JointPDF image info: " << metric->GetJointPDF() << std::endl;

  return EXIT_SUCCESS;
}
