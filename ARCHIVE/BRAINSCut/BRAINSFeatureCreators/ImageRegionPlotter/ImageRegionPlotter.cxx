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
// TODO Output Filename

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageRegionIterator.h"
#include "ImageRegionPlotterCLP.h"
#include "itkMinimumMaximumImageCalculator.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkHistogram.h"

#include "itkScalarImageToHistogramGenerator.h"
#include "itkHistogramToProbabilityImageFilter.h"
#include "itkMinimumMaximumImageCalculator.h"
#include "itkBinaryThresholdImageFilter.h"

#include "itkBRAINSROIAutoImageFilter.h"

#include "itkLabelStatisticsImageFilter.h"
#include "itkNumberToString.h"
#include <map>

#include <fstream>
#include <BRAINSCommonLib.h>

#define MAXIMUMLABELNUMBER 20
#define MAXIMUM_NUMBER_OF_INTENSITY 5000
/*
 * Author : Eun Young (Regina) Kim
 * Read in three types of images: Image, ROI mask, deformed probability Maps
 * Output will be...
 *    a text document with intensity values in it.
 */

/*
 * Spliting Filename to Directory and Filename itself
 * Ref: http://www.cplusplus.com/reference/string/string/find_last_of/
 */

/*
 *  get Path function
*/
std::string
GetPath( const std::string& filename)
{
  size_t found;

  found = filename.find_last_of("/");
  return filename.substr( 0, found );
}

// ///////////////////////////////////////////////////////////////////////////////////////
/*
 *  get Filename function
 */
std::string
GetFilename( const std::string& filename)
{
  size_t found;

  found = filename.find_last_of("/");
  return filename.substr( found + 1 );
}

// ///////////////////////////////////////////////////////////////////////////////////////
/*
 * Histogram Mapping to Image
 */
template <typename TImageType, typename THistogram>
typename TImageType::Pointer
MapHistogramToImage( typename TImageType::Pointer inputImage,
                     typename THistogram::Pointer histogram,
                     std::string mapFilename,
                     bool verbose)
{
  itk::NumberToString<double> doubleToString;

  /* ------------------------------------------------------------------------
   * Make Look Up Table
     ------------------------------------------------------------------------ */

  // Fine Maximum of Intensity

  using MaxCalculatorType = typename itk::MinimumMaximumImageCalculator<TImageType>;

  typename MaxCalculatorType::Pointer minmaxFilter = MaxCalculatorType::New();

  minmaxFilter->SetImage( inputImage );
  minmaxFilter->Compute();
  const int Image_MAXIMUM = static_cast<int>(ceil( minmaxFilter->GetMaximum() ) );
  const int Image_MINIMUM = static_cast<int>(ceil( minmaxFilter->GetMinimum() ) );

  std::cout << "Maximum Image Intensity :" << Image_MAXIMUM << std::endl;
  // Fixed Number of Maximum Number of Intensity
  if( (Image_MAXIMUM - Image_MINIMUM) > MAXIMUM_NUMBER_OF_INTENSITY )
    {
    itkGenericExceptionMacro( << "Range of Values are Too Large! "
                              << Image_MINIMUM << " " << Image_MAXIMUM);
    }

  // Compute Quantiles

  // - Initialize Intensity Look Up Table

  // TODO [lookUpTable] boolean variable used, produce look up table
  double lookUpTable[MAXIMUM_NUMBER_OF_INTENSITY] = {100.0};

  // - Generate Look Up Table

  std::ofstream mapOutFileStream;
  mapOutFileStream.open(mapFilename.c_str() );

  mapOutFileStream << "intensity, quantile "
                   << std::endl;

  int intensity = Image_MINIMUM;
  for( double quantile = 0.0;
       quantile < 1.0F;
       quantile = quantile + 0.01F )
    {
    double intensityAtQuantile = histogram->Quantile( 0, quantile );
    while( intensity < ceil( intensityAtQuantile ) )
      {
      lookUpTable[intensity] = quantile * 100;
      intensity++;
      }

    if( verbose )
      {
      std::cout << doubleToString(intensity) << " <--> "
                << doubleToString(quantile) << std::endl;
      }

    mapOutFileStream << doubleToString(intensity)
                     << "," << doubleToString(quantile)
                     << std::endl;
    }
  mapOutFileStream.close();

  /* ------------------------------------------------------------------------
   * - Generate Mapped ( Intensity->Quantile) Image
     ------------------------------------------------------------------------*/

  typename TImageType::Pointer outputImage = TImageType::New();
  outputImage->SetRegions( inputImage->GetLargestPossibleRegion() );
  outputImage->SetSpacing( inputImage->GetSpacing() );
  outputImage->SetOrigin( inputImage->GetOrigin() );
  outputImage->SetDirection( inputImage->GetDirection() );
  outputImage->Allocate();
  outputImage->FillBuffer( 0.0);

  using IteratorType = typename itk::ImageRegionIterator<TImageType>;

  IteratorType input_iterator( inputImage,
                               inputImage->GetLargestPossibleRegion() );

  IteratorType output_iterator( outputImage,
                                inputImage->GetLargestPossibleRegion() );
  for( input_iterator.GoToBegin(), output_iterator.GoToBegin();
       !input_iterator.IsAtEnd();
       ++input_iterator, ++output_iterator )
    {
    typename TImageType::PixelType currentIntensity = input_iterator.Get();
    output_iterator.Set( lookUpTable[(int)floor( currentIntensity )] );
    }

  return outputImage;
}

// ///////////////////////////////////////////////////////////////////////////////////////
/*
 * main function
 */
int
main(int argc, char *argv[])
{
  PARSE_ARGS;
  BRAINSRegisterAlternateIO();
  // define image with type of voxel
  using PixelType = double;
  constexpr int Dimension = 3;
  using ImageType = itk::Image<PixelType, Dimension>;

  // there has to be two input volumes and label volume
  if( (!inputLabelVolume.empty() ) && ( !inputVolume1.empty() ) &&
      (!inputVolume2.empty() ) )
    {
    // print volume names

    std::cout << "* InputImage Filename "  << inputVolume1  << std::endl
              << "* InputImage Filename "  << inputVolume2  << std::endl
              << "* Label Image Filename : " << inputLabelVolume << std::endl;

    // Image Reader for inputVolumes
    using ImageReaderType = itk::ImageFileReader<ImageType>;

    ImageReaderType::Pointer imageReader1 = ImageReaderType::New();
    imageReader1->SetFileName( inputVolume1 );

    imageReader1->Update();

    ImageReaderType::Pointer imageReader2 = ImageReaderType::New();
    imageReader2->SetFileName( inputVolume2 );

    imageReader2->Update();

    // Rescale input images between 0-4095

    using ProcessingPixelType = unsigned int;
    using ProcessingImageType = itk::Image<ProcessingPixelType, Dimension>;
    using RescalerType = itk::RescaleIntensityImageFilter<ImageType, ProcessingImageType>;

    constexpr unsigned int MIN = 0;
    constexpr unsigned int MAX = 4095;

    RescalerType::Pointer image1Rescaler = RescalerType::New();

    image1Rescaler->SetInput( imageReader1->GetOutput() );
    image1Rescaler->SetOutputMinimum( MIN );
    image1Rescaler->SetOutputMaximum( MAX );

    RescalerType::Pointer image2Rescaler = RescalerType::New();

    image2Rescaler->SetInput( imageReader2->GetOutput() );
    image2Rescaler->SetOutputMinimum( MIN );
    image2Rescaler->SetOutputMaximum( MAX );

    // read in input label image
    using LabelReaderType = itk::ImageFileReader<ProcessingImageType>;
    LabelReaderType::Pointer labelReader = LabelReaderType::New();

    labelReader->SetFileName( inputLabelVolume );

    labelReader->Update();

    // -  Get rid of Background
    // Get Brain Mask if given or call BRAINS ROI

    // - binary image type used for ROI

    using BinaryPixelType = unsigned char;
    using ROIImageType = itk::Image<BinaryPixelType, Dimension>;

    ROIImageType::ConstPointer roiVolume;

    // case 1: given inputBinaryROIVolume

    if( !inputBinaryROIVolume.empty() )
      {
      std::cout << "* Input Binary ROI Filename "  << inputBinaryROIVolume << std::endl;
      using ROIReaderType = itk::ImageFileReader<ROIImageType>;
      ROIReaderType::Pointer roiReader = ROIReaderType::New();

      roiReader->SetFileName( inputBinaryROIVolume );
      roiReader->Update();

      using ThresholdType = itk::BinaryThresholdImageFilter<ROIImageType, ROIImageType>;
      ThresholdType::Pointer thresholder = ThresholdType::New();

      thresholder->SetInput( roiReader->GetOutput() );
      thresholder->SetLowerThreshold( 1 );
      thresholder->SetOutsideValue(0);
      thresholder->SetInsideValue(1);
      thresholder->Update();

      roiVolume = thresholder->GetOutput();
      }
    else
      {
      // case 2: use ROIAUTO without any given binary image
      using ROIAutoType = itk::BRAINSROIAutoImageFilter<ProcessingImageType,
                                            ROIImageType>;

      ROIAutoType::Pointer image1_ROIFilter = ROIAutoType::New();

      image1_ROIFilter->SetInput( image1Rescaler->GetOutput() );

      image1_ROIFilter->Update();

      roiVolume = image1_ROIFilter->GetOutput();
      }

    // * Map Intensity Values to the Quantile

    ProcessingImageType::Pointer mapper1;
    ProcessingImageType::Pointer mapper2;

    if( !useIntensityForHistogram )
      {
      // typedef itk::Statistics::ScalarImagePortionToHistogramGenerator<
      // ProcessingImageType ,
      //                                                               MaskType>
      //                                                 histogramGeneratorType;
      using histogramGeneratorType = itk::LabelStatisticsImageFilter<ProcessingImageType,
                                              ROIImageType>;

      histogramGeneratorType::Pointer histogramGenerator1 = histogramGeneratorType::New();

      histogramGenerator1->SetInput( image1Rescaler->GetOutput() );
      histogramGenerator1->SetLabelInput( roiVolume );

      histogramGenerator1->SetHistogramParameters( numberOfHistogramBins,
                                                   MIN,
                                                   MAX );
      histogramGenerator1->UseHistogramsOn();
      histogramGenerator1->Update();

      histogramGeneratorType::Pointer histogramGenerator2 =
        histogramGeneratorType::New();

      histogramGenerator2->SetInput( image2Rescaler->GetOutput() );
      histogramGenerator2->SetLabelInput( roiVolume );
      histogramGenerator2->SetHistogramParameters( numberOfHistogramBins,
                                                   MIN,
                                                   MAX);
      histogramGenerator2->UseHistogramsOn();

      histogramGenerator2->Update();
      // * Generate Mapping Table from Intensity to Quantile
      ProcessingImageType::Pointer image1 = image1Rescaler->GetOutput();

      histogramGeneratorType::HistogramPointer histogram1 =
        histogramGenerator1->GetHistogram( 1 );

      mapper1 = MapHistogramToImage<ProcessingImageType,
                                    histogramGeneratorType::HistogramType>
          ( image1,
          histogram1,
          outputJointHistogramData + "map1.txt",
          verbose);

      ProcessingImageType::Pointer image2 = image2Rescaler->GetOutput();

      histogramGeneratorType::HistogramPointer histogram2 =
        histogramGenerator2->GetHistogram( 1 );

      mapper2 = MapHistogramToImage<ProcessingImageType,
                                    histogramGeneratorType::HistogramType>
          ( image2,
          histogram2,
          outputJointHistogramData + "map2.txt",
          verbose);
      }
    // * Iterator For Label Map

    // TODO [useIntensityForHistogram] used, produce intensity based histogram
    // as well
    using ConstIteratorType = itk::ImageRegionIterator<ProcessingImageType>;
    ConstIteratorType label_iterator( labelReader->GetOutput(),
                                      labelReader->GetOutput()->GetLargestPossibleRegion() );

    // Rescale Input Images
    using OutputImageType = itk::Image<unsigned int, 3>;
    using RescaleFilterType = itk::RescaleIntensityImageFilter<ProcessingImageType, OutputImageType>;

    RescaleFilterType::Pointer rescaler1 = RescaleFilterType::New();

    // Produce Joint Histogram of Intensity Values

    if( useIntensityForHistogram )
      {
      rescaler1->SetInput( image1Rescaler->GetOutput() );
      } // Produce Joint Histogram of Quantiles[ default ]
    else
      {
      rescaler1->SetInput( mapper1);
      }

    rescaler1->SetOutputMaximum(100);
    rescaler1->SetOutputMinimum(0);
    rescaler1->Update();

    RescaleFilterType::Pointer rescaler2 = RescaleFilterType::New();

    // Produce Joint Histogram of Intensity Values

    if( useIntensityForHistogram )
      {
      rescaler2->SetInput( image2Rescaler->GetOutput() );
      } // Produce Joint Histogram of Quantiles[ default ]
    else
      {
      rescaler2->SetInput( mapper2);
      }

    rescaler2->SetOutputMaximum(100);
    rescaler2->SetOutputMinimum(0);
    rescaler2->Update();

    // Start Iterator
    using LabelMapType = std::map<int, int>;
    LabelMapType labelIndex;
    int          new_key = 0;

    int HistogramArray[MAXIMUMLABELNUMBER][101][101] = { {{0}} };
    for( label_iterator.GoToBegin();
         !label_iterator.IsAtEnd();
         ++label_iterator )
      {
      if( label_iterator.Get() != 0 )
        {
        int temp_label =  label_iterator.Get();
        if( labelIndex.find( temp_label ) == labelIndex.end() )
          {
          if( verbose )
            {
            std::cout << "New Key Found : " << temp_label << std::endl;
            }
          labelIndex[temp_label] = new_key;
          new_key++;
          }

        int temp_image1_intensity =
          rescaler1->GetOutput()->GetPixel( label_iterator.GetIndex()  );
        int temp_image2_intensity =
          rescaler2->GetOutput()->GetPixel( label_iterator.GetIndex()  );

        HistogramArray[labelIndex[temp_label]][temp_image1_intensity][temp_image2_intensity]++;

        /*if( verbose ){
        std::cout << temp_label <<", "
                  << label_iterator.GetIndex() << " , "
                  << label_iterator.Get() << " , "
                  << temp_image1_intensity
                  <<", "
                  <<temp_image2_intensity
                  << std::endl;
        }*/
        }
      }

    // - open file stream and write to the file
    std::ofstream outputFileStream;
    outputFileStream.open( outputJointHistogramData.c_str() );

    const std::string IMAGE_FILENAME1 = GetFilename( inputVolume1);
    const std::string IMAGE_FILENAME2 = GetFilename( inputVolume2);
    // const std::string LABEL_MAP_FILENAME = GetFilename( inputLabelVolume );

    // - write header line 1: including image names

    outputFileStream << "[Image1]: " << inputVolume1 << std::endl
                     << "[Image2]: " << inputVolume2 << std::endl
                     << "[LabelMap]: " << inputLabelVolume
                     << std::endl;
    // - write header line 2: colume name

    std::string FrequencyName = "Frequence";
    if( useIntensityForHistogram )
      {
      FrequencyName += "OfIntensity";
      }
    else
      {
      FrequencyName += "OfQuantile";
      }

    outputFileStream << "label, " << IMAGE_FILENAME1
                     << ", " << IMAGE_FILENAME2
                     << ", " << FrequencyName
                     << std::endl;
    // - Iterate for each label type
    for( LabelMapType::iterator iter = labelIndex.begin();
         iter != labelIndex.end();
         ++iter )
      {
      // - Iterate for each bins
      for( int i = 0; i < 101; i++ )
        {
        for( int j = 0; j < 101; j++ )
          {
          // - Write text file

          outputFileStream <<  iter->first << ","
                           << i             << ","
                           << j             << ","
                           << HistogramArray[iter->second][i][j]
                           << std::endl;
          }
        }
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
