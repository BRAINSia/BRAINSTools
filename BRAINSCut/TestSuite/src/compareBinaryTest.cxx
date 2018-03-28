#include "itkImageFileReader.h"
#include "itkLabelOverlapMeasuresImageFilter.h"
#include "itkImage.h"
int main( int argc, char ** argv )
{
  int myExit = EXIT_SUCCESS;
  if( argc < 3 )
    {
    std::cerr << "Usage: " << std::endl;
    std::cerr << argv[0] << " inputBinaryFilename1 inputBinaryFilename2" << std::endl;
    myExit = EXIT_FAILURE;
    }


  typedef unsigned int PixelType;
  constexpr unsigned int Dimension = 3;
  typedef itk::Image< PixelType, Dimension >    LabelType;

  typedef itk::ImageFileReader< LabelType >  ReaderType;

  ReaderType::Pointer reader1 = ReaderType::New();
  ReaderType::Pointer reader2 = ReaderType::New();


  const char * inputBinaryFilename1  = argv[1];
  const char * inputBinaryFilename2  = argv[2];


  reader1->SetFileName( inputBinaryFilename1 );
  reader2->SetFileName( inputBinaryFilename2 );

  typedef itk::LabelOverlapMeasuresImageFilter< LabelType > OverlapFilterType;
  OverlapFilterType::Pointer filter = OverlapFilterType::New();

  filter->SetSourceImage( reader1->GetOutput() );
  filter->SetTargetImage( reader2->GetOutput() );

  try
    {
    filter->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    std::cerr << "ExceptionObject caught !" << std::endl;
    std::cerr << err << std::endl;
    myExit = EXIT_FAILURE;
    }

  if( filter->GetDiceCoefficient() < 0.97F )
  {
    std::cout<<"DSC = "<<  filter->GetDiceCoefficient() << std::endl;
    myExit = EXIT_FAILURE;
  }

  return myExit;
}
