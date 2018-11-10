#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkHessianToObjectnessMeasureImageFilter.h>
#include <itkMultiScaleHessianBasedMeasureImageFilter.h>
#include <itkRescaleIntensityImageFilter.h>
int main( int argc, char* argv[] )
{
  if( argc < 3 )
  {
    std::cerr << "Usage: "<< std::endl;
    std::cerr << argv[0];
    std::cerr << " <InputFileName> <OutputFileName>";
    std::cerr << " [SigmaMinimum] [SigmaMaximum] [NumberOfSigmaSteps]";
    std::cerr << " [alpha] [beta] [gamma] [structure] [bright]";
    std::cerr << std::endl;
    return EXIT_FAILURE;
  }
  const char * inputFileName = argv[1];
  const char * outputFileName = argv[2];
  double sigmaMinimum = 1.0;
  double sigmaMaximum = 10.0;
  unsigned int numberOfSigmaSteps = 10;
  double alpha = 0.1;
  double beta = 0.1;
  double gamma = 20.0;
  int structure = 2;
  bool bright = true;
  if( argc > 3 )
  {
    sigmaMinimum = atof( argv[3] );
  }
  if( argc > 4 )
  {
    sigmaMaximum = atof( argv[4] );
  }
  if( argc > 5 )
  {
    numberOfSigmaSteps = atof( argv[5] );
  }
  if( argc > 6 )
  {
    alpha = atof( argv[6] );
  }
  if( argc > 7 )
  {
    beta = atof(argv[7] );
  }
  if( argc > 8 )
  {
    gamma = atof( argv[8] );
  }
  if( argc > 9 )
  {
    structure = std::stoi( argv[9] );
  }
  if( argc >10)
  {
    bright = argv[10];
  }

  constexpr unsigned int Dimension = 3;
  using PixelType = float;
  using ImageType = itk::Image< PixelType, Dimension >;
  using HessianPixelType = itk::SymmetricSecondRankTensor< double, Dimension >;

  using HessianImageType = itk::Image< HessianPixelType, Dimension >;
  using ReaderType = itk::ImageFileReader< ImageType >;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( inputFileName );

  using ObjectnessFilterType = itk::HessianToObjectnessMeasureImageFilter< HessianImageType, ImageType >;
  ObjectnessFilterType::Pointer objectnessFilter = ObjectnessFilterType::New();
  objectnessFilter->SetBrightObject( bright );
  objectnessFilter->SetScaleObjectnessMeasure( false );
  objectnessFilter->SetAlpha( alpha );
  objectnessFilter->SetBeta( beta );
  objectnessFilter->SetGamma( gamma );
  objectnessFilter->SetObjectDimension(structure);

  using MultiScaleEnhancementFilterType = itk::MultiScaleHessianBasedMeasureImageFilter< ImageType, HessianImageType, ImageType >;
  MultiScaleEnhancementFilterType::Pointer multiScaleEnhancementFilter =
  MultiScaleEnhancementFilterType::New();
  multiScaleEnhancementFilter->SetInput( reader->GetOutput() );
  multiScaleEnhancementFilter->SetHessianToMeasureFilter( objectnessFilter );
  multiScaleEnhancementFilter->SetSigmaStepMethodToLogarithmic();
  multiScaleEnhancementFilter->SetSigmaMinimum( sigmaMinimum );
  multiScaleEnhancementFilter->SetSigmaMaximum( sigmaMaximum );
  multiScaleEnhancementFilter->SetNumberOfSigmaSteps( numberOfSigmaSteps );
  multiScaleEnhancementFilter->Update();

  using WriterType = itk::ImageFileWriter< ImageType >;
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( outputFileName );
  writer->SetInput( multiScaleEnhancementFilter->GetOutput() );
  std::cout<<objectnessFilter->GetAlpha()<<std::endl;
  std::cout<<objectnessFilter->GetBeta()<<std::endl;
  std::cout<<objectnessFilter->GetGamma()<<std::endl;


  try
  {
    writer->Update();
  }
  catch( itk::ExceptionObject & error )
  {
    std::cerr << "Error: " << error << std::endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
