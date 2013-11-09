// TODO Output Filename

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#include "itkOtsuHistogramMatchingImageFilter.h"
#include "itkHistogramMatchingImageFilter.h"

#include "HistogramMatchingFilterCLP.h"
#include "itkImageMaskSpatialObject.h"

#include <fstream>
#include <BRAINSCommonLib.h>

/*
 * Author : Eun Young (Regina) Kim modified by Hans J. Johnson
 *  Simple Histogram Matching Filter to see which values are good. :)
 *  This result is to be applied to the BRAINSCut procedure, including
 *  createVectors and applyModel for normalization purpose
 *
 */

/* -----------------------------------------------------------------------------
 * main function
 * ----------------------------------------------------------------------------- */

int
main(int argc, char *argv[])
{
  PARSE_ARGS;
  BRAINSRegisterAlternateIO();

  if( verbose )
    {
    std::cout << "- referenceVolume: " << referenceVolume << std::endl
              << "- inputVolume: " << inputVolume << std::endl
              << "- referenceBinaryVolume: " << referenceBinaryVolume << std::endl
              << "- inputBinaryVolume: " << inputBinaryVolume << std::endl
              << "- writeHistogram: " << writeHistogram << std::endl;
    }
  // define image with type of voxel
  typedef float PixelType;
  const unsigned int Dimension = 3;
  typedef itk::Image<PixelType, Dimension> ImageType;

  // Mask Reader type
  typedef itk::Image<unsigned char, Dimension> MaskVolumeType;
  typedef itk::ImageFileReader<MaskVolumeType> MaskReaderType;

  /*
   * Reader
   */
  typedef itk::ImageFileReader<ImageType> ImageReaderType;

  ImageReaderType::Pointer referenceImageReader = ImageReaderType::New();
  referenceImageReader->SetFileName( referenceVolume);

  ImageReaderType::Pointer inputImageReader = ImageReaderType::New();
  inputImageReader->SetFileName( inputVolume );

  /*
   * This part of program is for binary images of reference and input
   * image for histogram matching filter
   */

  typedef itk::ImageMaskSpatialObject<Dimension> MaskSpatialObjectType;
  MaskSpatialObjectType::Pointer referenceMaskSpatialObject = NULL;
  MaskSpatialObjectType::Pointer inputMaskSpatialObject = NULL;

  if( !referenceBinaryVolume.empty() && !inputBinaryVolume.empty() )
    {
    /*
     * Binary Image Read
     */
    typedef unsigned char                          BinaryPixelType;
    typedef itk::Image<BinaryPixelType, Dimension> BinaryImageType;

    typedef itk::ImageFileReader<BinaryImageType> BinaryImageReaderType;

    BinaryImageReaderType::Pointer referenceBinaryImageReader =
      BinaryImageReaderType::New();
    referenceBinaryImageReader->SetFileName( referenceBinaryVolume );
    referenceBinaryImageReader->Update();

    BinaryImageReaderType::Pointer inputBinaryImageReader =
      BinaryImageReaderType::New();
    inputBinaryImageReader->SetFileName( inputBinaryVolume );
    inputBinaryImageReader->Update();

    /*
     * Binary image to Spatial Object
     */
    if( verbose )
      {
      std::cout << " using binary masks in processing histograms " << std::endl;
      }
    referenceMaskSpatialObject = MaskSpatialObjectType::New();
    inputMaskSpatialObject = MaskSpatialObjectType::New();

    referenceMaskSpatialObject->SetImage( referenceBinaryImageReader->GetOutput() );
    inputMaskSpatialObject->SetImage( inputBinaryImageReader->GetOutput() );

    referenceMaskSpatialObject->ComputeObjectToWorldTransform();
    inputMaskSpatialObject->ComputeObjectToWorldTransform();

    if( referenceMaskSpatialObject.IsNotNull() )
      {
      std::cout << "ref is not null" << std::endl;
      }
    else
      {
      std::cout << "ref is null" << std::endl;
      }
    }

  /*
   *  Histogram Matching Filter
   */

  std::string histogramDataFilename = writeHistogram + ".dat";
  // Define Writer Here
  typedef itk::ImageFileWriter<ImageType> ImageWriterType;
  ImageWriterType::Pointer imageWriter = ImageWriterType::New();
  imageWriter->SetFileName( outputVolume );

  /*
   * Default Behavior and Otsu Behavior
   *  itk::HistogramMatchingImageFilter
   */

  if( histogramAlgorithm == "quantileHistogramMatch" )
    {
    std::cout << "Using: " << histogramAlgorithm << " algorithm." << std::endl;
    std::cout << "This is under the construction. " << std::endl;
    return EXIT_FAILURE;
    // Find quantiles and linearly scale image.
    }
  if( histogramAlgorithm == "simpleITKHistogramMatch" )
    {
    std::cout << "Using: " << histogramAlgorithm << " algorithm." << std::endl;
    typedef itk::HistogramMatchingImageFilter<ImageType,
                                              ImageType> SimpleHistogramMatchingType;
    SimpleHistogramMatchingType::Pointer SHFilter = SimpleHistogramMatchingType::New();

    SHFilter->SetReferenceImage( referenceImageReader->GetOutput() );
    SHFilter->SetInput( inputImageReader->GetOutput() );
    SHFilter->SetNumberOfHistogramLevels( numberOfHistogramBins );
    SHFilter->SetNumberOfMatchPoints( numberOfMatchPoints );
    SHFilter->ThresholdAtMeanIntensityOn();

    /*  Writer  */
    imageWriter->SetInput( SHFilter->GetOutput() );
    imageWriter->Update();
    // ---------------------------------------------------------------------------
    // //
    // To write out histogram data for plotting with R
    // - This part should write
    // --   Histogram data for reference, input, and output volume
    // --   Shell script to run R script.
    // --   The R scrip is included in the svn repository
    // "HistogramMatchingFilter.R"
    // ---------------------------------------------------------------------------
    // //

    typedef SimpleHistogramMatchingType::HistogramType HistogramType;

    HistogramType::ConstPointer srcHG = SHFilter->GetSourceHistogram();
    HistogramType::ConstPointer refHG = SHFilter->GetReferenceHistogram();
    HistogramType::ConstPointer outHG = SHFilter->GetOutputHistogram();

    HistogramType::ConstIterator srcIt = srcHG->Begin();
    HistogramType::ConstIterator refIt = refHG->Begin();
    HistogramType::ConstIterator outIt = outHG->Begin();

    HistogramType::ConstIterator endIt = srcHG->End();

    std::ofstream histFileStream;

    histFileStream.open( histogramDataFilename.c_str() );

    histFileStream << " bin, src_frequency, ref_frequency, adjusted_frequency"
                   << std::endl;

    unsigned int binNumber = 0;

    while( srcIt != endIt )
      {
      histFileStream << binNumber << ", "
                     << srcIt.GetFrequency() << ", "
                     << refIt.GetFrequency() << ", "
                     << outIt.GetFrequency()
                     << std::endl;
      ++srcIt;
      ++refIt;
      ++outIt;
      binNumber++;
      }

    histFileStream.close();
    }
  else if( histogramAlgorithm == "OtsuHistogramMatching" )
    {
    std::cout << "Using: " << histogramAlgorithm << " algorithm." << std::endl;
    typedef itk::OtsuHistogramMatchingImageFilter<ImageType,
                                                  ImageType> OtsuHistogramMatchingType;
    OtsuHistogramMatchingType::Pointer OHFilter = OtsuHistogramMatchingType::New();

    OHFilter->SetReferenceImage( referenceImageReader->GetOutput() );
    OHFilter->SetInput( inputImageReader->GetOutput() );

    if( referenceMaskSpatialObject.IsNotNull() )
      {
      std::cout << "Setting Reference Mask" << std::endl;
      OHFilter->SetReferenceMask( referenceMaskSpatialObject );
      }
    if( inputMaskSpatialObject.IsNotNull() )
      {
      std::cout << "Setting Input Mask" << std::endl;
      OHFilter->SetSourceMask(inputMaskSpatialObject);
      }

    OHFilter->SetNumberOfHistogramLevels( numberOfHistogramBins );
    OHFilter->SetNumberOfMatchPoints( numberOfMatchPoints );

    /* Writer */
    imageWriter->SetInput( OHFilter->GetOutput() );
    imageWriter->Update();
    // ---------------------------------------------------------------------------
    // //
    // To write out histogram data for plotting with R
    // - This part should write
    // --   Histogram data for reference, input, and output volume
    // --   Shell script to run R script.
    // --   The R scrip is included in the svn repository
    // "HistogramMatchingFilter.R"
    // ---------------------------------------------------------------------------
    // //

    typedef  OtsuHistogramMatchingType::HistogramType HistogramType;
    HistogramType::ConstPointer srcHG = OHFilter->GetSourceHistogram();
    HistogramType::ConstPointer refHG = OHFilter->GetReferenceHistogram();
    HistogramType::ConstPointer outHG = OHFilter->GetOutputHistogram();

    HistogramType::ConstIterator srcIt = srcHG->Begin();
    HistogramType::ConstIterator refIt = refHG->Begin();
    HistogramType::ConstIterator outIt = outHG->Begin();
    HistogramType::ConstIterator endIt = srcHG->End();

    std::ofstream histFileStream;
    histFileStream.open( histogramDataFilename.c_str() );

    histFileStream << " bin, src_frequency, ref_frequency, adjusted_frequency"
                   << std::endl;

    unsigned int binNumber = 0;

    while( srcIt != endIt )
      {
      histFileStream << binNumber << ", "
                     << srcIt.GetFrequency() << ", "
                     << refIt.GetFrequency() << ", "
                     << outIt.GetFrequency()
                     << std::endl;
      ++srcIt;
      ++refIt;
      ++outIt;
      binNumber++;
      }

    histFileStream.close();
    }
  else
    {
    std::cout << "Unsupported Histogram Algorithm Option "
              << histogramAlgorithm
              << ". Please check available option with  -help"
              << std::endl;
    return EXIT_FAILURE;
    }

  /*
   * write shell script to run R script for plotting
   */

  std::ofstream histgramPlotterStream;
  std::string   histogramPloterFilename = writeHistogram + ".sh";
  histgramPlotterStream.open( histogramPloterFilename.c_str() );

  std::string histogramFigureFilename = writeHistogram + ".png";
  std::string histogramRScript = "HistogramMatchingFilter.R";

  histgramPlotterStream << " R --slave --args "
                        << histogramDataFilename << " "
                        << histogramFigureFilename << " "
                        << " < "
                        << histogramRScript
                        << std::endl;
  return EXIT_SUCCESS;
}
