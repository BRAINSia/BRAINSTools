// TODO Output Filename

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionIterator.h"
#include "JointHistogramCLP.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkHistogram.h"

#include "itkScalarImageToHistogramGenerator.h"
#include "itkHistogramToProbabilityImageFilter.h"
#include "itkMinimumMaximumImageCalculator.h"
#include "itkInvertIntensityImageFilter.h"

#include "itkLabelStatisticsImageFilter.h"
#include <map>

#include <fstream>

#define HISTOGRAMSIZE 50

/*
 * Author : Eun Young (Regina) Kim
 */

// ///////////////////////////////////////////////////////////////////////////////////////
/*
 * main function
 */
int
main(int argc, char *argv[])
{
  PARSE_ARGS;
  // define image with type of voxel
  typedef double PixelType;
  const int Dimension = 3;
  typedef itk::Image<PixelType, Dimension> InputImageType;

  // there has to be two input volumes and label volume
  if( ( !inputVolume1.empty() ) && (!inputVolume2.empty() ) )
    {
    // print volume names

    std::cout << "* InputImage Filename "  << inputVolume1  << std::endl
              << "* InputImage Filename "  << inputVolume2  << std::endl;

    // Image Reader for inputVolumes
    typedef itk::ImageFileReader<InputImageType> ImageReaderType;

    ImageReaderType::Pointer imageReader1 = ImageReaderType::New();
    imageReader1->SetFileName( inputVolume1 );

    ImageReaderType::Pointer imageReader2 = ImageReaderType::New();
    imageReader2->SetFileName( inputVolume2 );
    try
      {
      imageReader1->Update();
      imageReader2->Update();
      }
    catch( itk::ExceptionObject & e )
      {
      std::cerr << "Exception in Resampling." << std::endl;
      std::cerr << e.GetDescription() << std::endl;
      std::cerr << e.GetLocation() << std::endl;
      return EXIT_FAILURE;
      }

    // Rescale Input Images
    typedef itk::RescaleIntensityImageFilter<InputImageType, InputImageType>
      RescaleFilterType;

    RescaleFilterType::Pointer rescaler1 = RescaleFilterType::New();

    rescaler1->SetInput( imageReader1->GetOutput() );

    rescaler1->SetOutputMaximum( HISTOGRAMSIZE );
    rescaler1->SetOutputMinimum(0);

    RescaleFilterType::Pointer rescaler2 = RescaleFilterType::New();

    rescaler2->SetInput( imageReader2->GetOutput() );
    rescaler2->SetOutputMaximum(HISTOGRAMSIZE);
    rescaler2->SetOutputMinimum(0);
    try
      {
      rescaler1->Update();
      rescaler2->Update();
      }
    catch( itk::ExceptionObject & e )
      {
      std::cerr << "Exception in Resampling." << std::endl;
      std::cerr << e.GetDescription() << std::endl;
      std::cerr << e.GetLocation() << std::endl;
      return EXIT_FAILURE;
      }

    // Start Iterator
    unsigned int HistogramArray[HISTOGRAMSIZE][HISTOGRAMSIZE] =  {{0}};
    // * Iterator For Image

    typedef itk::ImageRegionIterator<InputImageType> ConstIteratorType;
    ConstIteratorType it( imageReader1->GetOutput(),
                          imageReader1->GetOutput()->GetLargestPossibleRegion() );
    for( it.GoToBegin();
         !it.IsAtEnd();
         ++it )
      {
      if( it.Get() != 0 )
        {
        int temp_image1_intensity =
          rescaler1->GetOutput()->GetPixel( it.GetIndex()  );
        int temp_image2_intensity =
          rescaler2->GetOutput()->GetPixel( it.GetIndex()  );

        HistogramArray[temp_image1_intensity][temp_image2_intensity]++;

        if( verbose )
          {
          std::cout << " ADD to the Bin (" << temp_image1_intensity
                    << " , " << temp_image2_intensity << " ) "
                    << " = " << HistogramArray[temp_image1_intensity][temp_image2_intensity]
                    << std::endl;
          }
        }
      }
    // - open file stream and write to the file
    std::ofstream outputFileStream;
    std::string   outputHistogramData = outputJointHistogramImage + ".txt";
    outputFileStream.open( outputHistogramData.c_str() );

    // - write header line 1: including image names

    outputFileStream << "[Image1]: " << inputVolume1 << std::endl
                     << "[Image2]: " << inputVolume2 << std::endl
                     << std::endl;
    // - write header line 2: colume name

    std::string FrequencyName = "Frequence";
    FrequencyName += "OfIntensity";

    outputFileStream << "label, " << inputVolume1
                     << ", " << inputVolume2
                     << ", " << FrequencyName
                     << std::endl;
    // - Iterate for each bins
    for( int i = 0; i < HISTOGRAMSIZE; i++ )
      {
      for( int j = 0; j < HISTOGRAMSIZE; j++ )
        {
        // - Write text file

        outputFileStream
          << i             << ","
          << j             << ","
          << HistogramArray[i][j]
          << std::endl;
        }
      }
    // - Write Histogram Image
    typedef itk::Image<unsigned int, 2> HistogramImageType;
    HistogramImageType::Pointer histogramImg = HistogramImageType::New();

    HistogramImageType::IndexType start;
    start[0] = 0; start[1] = 0;

    HistogramImageType::SizeType size;
    size[0] = HISTOGRAMSIZE; size[1] = HISTOGRAMSIZE;

    HistogramImageType::RegionType region;
    region.SetSize( size);
    region.SetIndex( start );

    histogramImg->SetRegions( region );

    HistogramImageType::SpacingType space;
    space[0] = 1; space[1] = 1;
    histogramImg->SetSpacing( space );

    histogramImg->Allocate();

    // * Iterator For Image

    typedef itk::ImageRegionIterator<HistogramImageType> HistIteratorType;
    HistIteratorType hit( histogramImg,
                          histogramImg->GetLargestPossibleRegion() );
    for( hit.GoToBegin(); !hit.IsAtEnd();     ++hit )
      {
      HistogramImageType::IndexType currentIdx = hit.GetIndex();
      hit.Set( HistogramArray[currentIdx[0]][currentIdx[1]] );
      }

    // Histogram Image Rescale to 0-255

    typedef itk::Image<unsigned char, 2> HistogramWritingType;
    typedef itk::RescaleIntensityImageFilter<HistogramImageType,
                                             HistogramWritingType>
      HistogramRescaleFilterType;

    HistogramRescaleFilterType::Pointer histogramRescaler
      = HistogramRescaleFilterType::New();

    histogramRescaler->SetInput( histogramImg );

    histogramRescaler->SetOutputMaximum(255);
    histogramRescaler->SetOutputMinimum(0);

    // Invert Intensity so that background is white
    typedef itk::InvertIntensityImageFilter<HistogramWritingType,
                                            HistogramWritingType>
      InvertIntensityFilterType;
    InvertIntensityFilterType::Pointer invertFilter =
      InvertIntensityFilterType::New();

    invertFilter->SetInput( histogramRescaler->GetOutput() );
    invertFilter->SetMaximum(255);
    invertFilter->Update();

    //  Histogram Writer
    typedef itk::ImageFileWriter<HistogramWritingType> HistogramWriter;

    HistogramWriter::Pointer histogramWriter  = HistogramWriter::New();

    histogramWriter->SetFileName( outputJointHistogramImage );
    histogramWriter->SetInput( invertFilter->GetOutput() );

    try
      {
      histogramWriter->Update();
      }
    catch( itk::ExceptionObject & e )
      {
      std::cerr << "Exception in Resampling." << std::endl;
      std::cerr << e.GetDescription() << std::endl;
      std::cerr << e.GetLocation() << std::endl;
      return EXIT_FAILURE;
      }
    }
  else
    {
    std::cout << " Wrong Argument! " << std::endl
              << " inputVolume1, inputVolume2, and inputLabelVolume are necessary! "
              << std::endl;
    return EXIT_FAILURE;
    }
  return 0;
}
