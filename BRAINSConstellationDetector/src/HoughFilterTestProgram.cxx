// Inspecting the performance of the HoughTransformRadialVotingImageFilter regarding the number of threads

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkHoughTransformRadialVotingImageFilter.h"

int main( int argc, char * argv[] )
{
  const unsigned int LocalImageDimension = 3;

  typedef short                                            InputPixelType;
  typedef double                                           OutputPixelType;
  typedef itk::Image<InputPixelType, LocalImageDimension>  InputImageType;
  typedef itk::Image<OutputPixelType, LocalImageDimension> OutputImageType;

  typedef itk::HoughTransformRadialVotingImageFilter<InputImageType, OutputImageType> HoughFilterType;
  typedef HoughFilterType::SpheresListType                                            SpheresListType;
  typedef HoughFilterType::Pointer                                                    HoughFilterPointer;
  HoughFilterPointer houghFilter = HoughFilterType::New();

// Reader type
  typedef itk::ImageFileReader<InputImageType> ReaderType;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[1] );
// reader->Update();

  houghFilter->SetInput( reader->GetOutput() );

  houghFilter->SetNumberOfSpheres( 2 );
  houghFilter->SetMinimumRadius( 11. );
  houghFilter->SetMaximumRadius( 13. );
  houghFilter->SetSigmaGradient( 1. );
  houghFilter->SetVariance( 1. );
  houghFilter->SetSphereRadiusRatio( 1. );
  houghFilter->SetVotingRadiusRatio( .5 );
  houghFilter->SetThreshold( 10. );
  houghFilter->SetOutputThreshold( .8 );
  houghFilter->SetGradientThreshold( 0. );
  houghFilter->SetNbOfThreads( 64 );
  houghFilter->SetSamplingRatio( .2 );
  houghFilter->SetHoughEyeDetectorMode( 1 );

  try
    {
    houghFilter->Update();
    }
  catch( itk::ExceptionObject & excep )
    {
    std::cerr << "Failed houghFilter " << std::endl;
    std::cerr << excep << std::endl;
    }
  catch( ... )
    {
    std::cout << "Failed on houghFilter exception occured" << std::endl;
    }

/*
this->m_AccumulatorImage = houghFilter->GetOutput();

// Write debug accumulator image
typename WriterType::Pointer writer = WriterType::New();
writer->SetFileName( this->m_ResultsDir + "/HoughEyeAccumulator.nii.gz" );
writer->SetInput( this->m_AccumulatorImage );
writer->SetUseCompression( true );

*/

// writer type
  typedef  itk::ImageFileWriter<OutputImageType> WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName(argv[2]);
  writer->SetInput( houghFilter->GetOutput() );
  try
    {
    writer->Update();
    }
  catch( itk::ExceptionObject & excep )
    {
    std::cerr << "Cannot write the Accumulator image!" << std::endl;
    std::cerr << excep << std::endl;
    }

  const SpheresListType spheres = houghFilter->GetSpheres();
  if( spheres.size() < 2 )
    {
    std::cerr << "Error: The number of detected spheres is less than 2!" << std::endl;
    // std::cerr << "The program will continue to run for generating some debug output for GUI corrector." << std::endl;
    // std::cerr << "Make sure you set the debug level > 4." << std::endl;
    // this->m_Failure = true;
    // return -1;
    }
  else
    {
    std::cout << "It works! The number of detected spheres is not less than 2!" << std::endl;
    }

  return EXIT_SUCCESS;
}
