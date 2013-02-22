#include "BRAINSCommonLib.h"

#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionIterator.h"
#include "itkTestingExtractSliceImageFilter.h"
#include "itkCastImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkTileImageFilter.h"
#include "itkFlipImageFilter.h"
#include "itkLabelOverlayImageFilter.h"
#include "itkComposeImageFilter.h"
#include "itkRGBPixel.h"

#include "BRAINSSnapShotWriterCLP.h"

/*
 * extracting slice numbers in index
 */

typedef std::vector<unsigned int> ExtractIndexType;

typedef std::vector<int>   IndexType;
typedef std::vector<int>   PercentIndexType;
typedef std::vector<float> PhysicalPointIndexType;

template <class TImageType>
ExtractIndexType  GetSliceIndexToExtract(
  typename TImageType::Pointer   referenceImage,
  std::vector<int>      planes,
  IndexType                      inputSliceToExtractInIndex,
  PercentIndexType               inputSliceToExtractInPercent,
  PhysicalPointIndexType          inputSliceToExtractInPhysicalPoint)
{
  if( inputSliceToExtractInIndex.empty() &&
      inputSliceToExtractInPercent.empty() &&
      inputSliceToExtractInPhysicalPoint.empty() )
    {
    std::cout << "ERROR:: one of input index has to be entered "
              << std::endl;
    exit(EXIT_FAILURE);
    }

  ExtractIndexType sliceIndexToExtract;
  if( !inputSliceToExtractInIndex.empty() )
    {
    for( unsigned int i = 0; i < inputSliceToExtractInIndex.size(); i++ )
      {
      sliceIndexToExtract.push_back( inputSliceToExtractInIndex[i] );
      }
    }
  else if( !inputSliceToExtractInPhysicalPoint.empty() )
    {
    for( unsigned int i = 0; i < inputSliceToExtractInPhysicalPoint.size(); i++ )
      {
      typename TImageType::PointType physicalPoints;
      typename TImageType::IndexType dummyIndex;
      for( unsigned int p = 0; p < physicalPoints.Size(); p++ )
        {
        // fill the same value
        physicalPoints[p] = inputSliceToExtractInPhysicalPoint[i];
        }
      referenceImage->TransformPhysicalPointToIndex( physicalPoints,
                                                     dummyIndex );

      std::cout << inputSliceToExtractInPhysicalPoint[i]
                << "-->"
                << dummyIndex[planes[i]]
                << std::endl;
      sliceIndexToExtract.push_back( dummyIndex[planes[i]] );
      }
    }
  else if( !inputSliceToExtractInPercent.empty() )
    {
    for( unsigned int i = 0; i < inputSliceToExtractInPercent.size(); i++ )
      {
      if( inputSliceToExtractInPercent[i] < 0.0F ||
          inputSliceToExtractInPercent[i] > 100.0F )
        {
        std::cout << "ERROR: Percent has to be between 0 and 100 "
                  << std::endl;
        exit( EXIT_FAILURE );
        }
      unsigned int size = (referenceImage->GetLargestPossibleRegion() ).GetSize()[planes[i]];
      unsigned int index =
        ( (float)inputSliceToExtractInPercent[i] / 100.0F ) * size;

      std::cout << inputSliceToExtractInPercent[i]
                << "-->"
                << index
                << std::endl;
      sliceIndexToExtract.push_back( index );
      }
    }

  return sliceIndexToExtract;
}

/*
 * change orientation
 */
template <class TImageType>
// input parameter type
typename TImageType::Pointer ChangeOrientOfImage( typename TImageType::Pointer imageVolume,
                                                  itk::FixedArray<bool, 3> flipAxes )
{
  typedef itk::FlipImageFilter<TImageType> FlipImageFilterType;

  typename FlipImageFilterType::Pointer flipFilter =
    FlipImageFilterType::New();

  flipFilter->SetInput( imageVolume );
  flipFilter->SetFlipAxes( flipAxes );
  try
    {
    flipFilter->Update();
    }
  catch( ... )
    {
    std::cout << "ERROR: Fail to flip the image "
              << std::endl;
    }

  return flipFilter->GetOutput();
}

/*
 * template reading function
 */
template <class TStringVectorType, // input parameter type
          class TReaderType,       // reader type
          class TImageVectorType>
// return type
TImageVectorType ReadImageVolumes( TStringVectorType filenameVector )
{
  typedef typename TReaderType::Pointer                  ReaderPointer;
  typedef typename TReaderType::OutputImageType::Pointer OutputImagePointerType;

  TImageVectorType imageVector;
  for( unsigned int i = 0; i < filenameVector.size(); i++ )
    {
    std::cout << "Reading image " << i + 1 << ": " << filenameVector[i] << "...\n";

    ReaderPointer reader = TReaderType::New();
    reader->SetFileName( filenameVector[i].c_str() );

    try
      {
      reader->Update();
      }
    catch( ... )
      {
      std::cout << "ERROR:  Could not read image " << filenameVector[i] << "." << std::endl;
      exit(EXIT_FAILURE);
      }

    OutputImagePointerType image = reader->GetOutput();

    itk::FixedArray<bool, 3> flipAxes;
    flipAxes[0] = 0;
    flipAxes[1] = 0;
    flipAxes[2] = 1;

    OutputImagePointerType orientedImage =
      ChangeOrientOfImage<typename TReaderType::OutputImageType>( image, flipAxes );

    imageVector.push_back( orientedImage );
    }

  return TImageVectorType( imageVector );
}

/*
 * extract slices
 */
template <class TInputImageType,
          class TOutputImageType>
typename TOutputImageType::Pointer
ExtractSlice( typename TInputImageType::Pointer inputImage,
              int plane,
              int sliceNumber)
{
  if( plane < 0 || plane > 3 )
    {
    std::cout << "ERROR: Extracting plane should be between 0 and 2(0,1,or 2)" << std::endl;
    exit(EXIT_FAILURE);
    }
  /* extract 2D plain */
  typedef itk::Testing::ExtractSliceImageFilter<TInputImageType,
                                                TOutputImageType> ExtractVolumeFilterType;

  typename ExtractVolumeFilterType::Pointer extractVolumeFilter = ExtractVolumeFilterType::New();

  typename TInputImageType::RegionType region = inputImage->GetLargestPossibleRegion();

  typename TInputImageType::SizeType size = region.GetSize();
  size[plane] = 0;

  typename TInputImageType::IndexType start = region.GetIndex();
  start[plane] = sliceNumber;

  typename TInputImageType::RegionType outputRegion;
  outputRegion.SetSize( size );
  outputRegion.SetIndex( start );

  extractVolumeFilter->SetExtractionRegion( outputRegion );
  extractVolumeFilter->SetInput( inputImage );
  extractVolumeFilter->SetDirectionCollapseToGuess();
  extractVolumeFilter->Update();

  typename TOutputImageType::Pointer outputImage = extractVolumeFilter->GetOutput();
  return outputImage;
}

/* scaling between 0-255 */
template <class TInputImage, class TOutputImage>
typename TOutputImage::Pointer
Rescale( const typename TInputImage::Pointer inputImage,
         const int min,
         const int max)
{
  typedef itk::RescaleIntensityImageFilter<TInputImage,
                                           TInputImage> RescaleFilterType;

  typename RescaleFilterType::Pointer rescaler = RescaleFilterType::New();

  rescaler->SetInput( inputImage );
  rescaler->SetOutputMinimum( min );
  rescaler->SetOutputMaximum( max );

  typedef typename itk::CastImageFilter<TInputImage,
                                        TOutputImage> CastingFilterType;

  typename CastingFilterType::Pointer caster = CastingFilterType::New();
  caster->SetInput( rescaler->GetOutput() );
  caster->Update();

  typename TOutputImage::Pointer outputImage = caster->GetOutput();

  return outputImage;
}

/*
 * main
 */

int
main(int argc, char * *argv)
{
  PARSE_ARGS;

  if( inputVolumes.empty() )
    {
    std::cout << "Input image volume is required "
              << std::endl;
    exit(EXIT_FAILURE);
    }
  if( inputPlaneDirection.size() == 0 )
    {
    std::cout << "Input Plane Direction is required "
              << std::endl;
    exit(EXIT_FAILURE);
    }
  if( inputSliceToExtractInIndex.size() == 0 &&
      inputSliceToExtractInPercent.size() == 0 &&
      inputSliceToExtractInPhysicalPoint.size() )
    {
    std::cout << "At least one of input Slice to Extract has to be specified."
              << std::endl;
    exit(EXIT_FAILURE);
    }

  if( inputPlaneDirection.size() != inputSliceToExtractInIndex.size() &&
      inputPlaneDirection.size() != inputSliceToExtractInPercent.size() &&
      inputPlaneDirection.size() != inputSliceToExtractInPhysicalPoint.size() )
    {
    std::cout << "Number of input slice number should be equal input plane direction."
              << std::endl;
    exit(EXIT_FAILURE);
    }

  const unsigned int numberOfImgs = inputVolumes.size();

  /* type definition */
  typedef itk::Image<double, 3>        Image3DVolumeType;
  typedef itk::Image<double, 2>        Image2DVolumeType;
  typedef itk::Image<unsigned char, 3> Image3DBinaryType;

  typedef std::vector<std::string>                ImageFilenameVectorType;
  typedef std::vector<Image3DVolumeType::Pointer> Image3DVolumeVectorType;
  typedef std::vector<Image3DBinaryType::Pointer> Image3DBinaryVectorType;

  typedef itk::ImageFileReader<Image3DVolumeType> Image3DVolumeReaderType;
  typedef Image3DVolumeReaderType::Pointer        Image3DVolumeReaderPointer;

  typedef itk::ImageFileReader<Image3DBinaryType> Image3DBinaryReaderType;
  typedef Image3DBinaryReaderType::Pointer        Image3DBinaryReaderPointer;

  typedef itk::Image<unsigned char, 2> OutputGreyImageType;

  typedef itk::RGBPixel<unsigned char> RGBPixelType;
  typedef itk::Image<RGBPixelType, 2>  OutputRGBImageType;

  /* read in image volumes */
  Image3DVolumeVectorType image3DVolumes = ReadImageVolumes<ImageFilenameVectorType,
                                                            Image3DVolumeReaderType,
                                                            Image3DVolumeVectorType>
      ( inputVolumes );

  /* read in binary volumes */
  Image3DBinaryVectorType image3DBinaries = ReadImageVolumes<ImageFilenameVectorType,
                                                             Image3DBinaryReaderType,
                                                             Image3DBinaryVectorType>
      ( inputBinaryVolumes );

  ExtractIndexType extractingSlices =
    GetSliceIndexToExtract<Image3DVolumeType>( image3DVolumes[0],
                                               inputPlaneDirection,
                                               inputSliceToExtractInIndex,
                                               inputSliceToExtractInPercent,
                                               inputSliceToExtractInPhysicalPoint);

  /* combine binary images */
  Image3DVolumeType::Pointer labelMap = Image3DVolumeType::New();
  if( !image3DBinaries.empty() )
    {
    labelMap->CopyInformation( image3DBinaries[0] );
    labelMap->SetRegions( image3DVolumes[0]->GetLargestPossibleRegion() );
    labelMap->Allocate();

    itk::ImageRegionIterator<Image3DVolumeType> binaryIterator(
      labelMap,
      labelMap->GetLargestPossibleRegion() );

    while( !binaryIterator.IsAtEnd() )
      {
      Image3DBinaryType::IndexType index = binaryIterator.GetIndex();
      /** itereate one image */
      binaryIterator.Set( 0 );
      for( unsigned int i = 0; i < image3DBinaries.size(); i++ )
        {
        if( image3DBinaries[i]->GetPixel( index ) > 0 )
          {
          binaryIterator.Set( i + 1 ); // label color zero is grey
          }
        }
      /** add each binary values */
      ++binaryIterator;
      }
    }
  /* compose color image */
  typedef itk::LabelOverlayImageFilter<OutputGreyImageType,
                                       OutputGreyImageType,
                                       OutputRGBImageType> LabelOverlayFilter;

  typedef itk::ComposeImageFilter<OutputGreyImageType,
                                  OutputRGBImageType> RGBComposeFilter;

  typedef std::vector<OutputRGBImageType::Pointer> OutputRGBImageVectorType;

  OutputRGBImageVectorType rgbSlices;
  for( unsigned int plane = 0; plane < inputPlaneDirection.size(); plane++ )
    {
    for( unsigned int i = 0; i < numberOfImgs; i++ )
      {
      /** get slicer */
      Image3DVolumeType::Pointer current3DImage = image3DVolumes[i];
      Image2DVolumeType::Pointer imageSlice =
        ExtractSlice<Image3DVolumeType, Image2DVolumeType>( current3DImage,
                                                            inputPlaneDirection[plane],
                                                            extractingSlices[plane] );

      OutputGreyImageType::Pointer greyScaleSlice =
        Rescale<Image2DVolumeType, OutputGreyImageType>( imageSlice, 0, 255 );

      /** binaries */
      OutputGreyImageType::Pointer labelSlice;
      if( !image3DBinaries.empty() )
        {
        /** rgb creator */
        LabelOverlayFilter::Pointer rgbComposer = LabelOverlayFilter::New();

        labelSlice =
          ExtractSlice<Image3DVolumeType, OutputGreyImageType>( labelMap,
                                                                inputPlaneDirection[plane],
                                                                extractingSlices[plane] );
        rgbComposer->SetLabelImage( labelSlice);
        rgbComposer->SetInput( greyScaleSlice );
        rgbComposer->SetOpacity(.5F);
        try
          {
          rgbComposer->Update();
          }
        catch( itk::ExceptionObject& e  )
          {
          std::cout << "ERROR:  Could not update image." << std::endl;
          std::cout << "ERROR:  " << e.what() << std::endl;
          exit(EXIT_FAILURE);
          }

        rgbSlices.push_back( rgbComposer->GetOutput() );
        }
      else /** ----------------------------------------------------- */
        {
        RGBComposeFilter::Pointer rgbComposer = RGBComposeFilter::New();

        rgbComposer->SetInput1( greyScaleSlice );
        rgbComposer->SetInput2( greyScaleSlice );
        rgbComposer->SetInput3( greyScaleSlice );

        try
          {
          rgbComposer->Update();
          }
        catch( itk::ExceptionObject& e  )
          {
          std::cout << "ERROR:  Could not update image." << std::endl;
          std::cout << "ERROR:  " << e.what() << std::endl;
          exit(EXIT_FAILURE);
          }
        rgbSlices.push_back( rgbComposer->GetOutput() );
        }
      }
    }

  /* tile the images */
  typedef itk::TileImageFilter<OutputRGBImageType, OutputRGBImageType> TileFilterType;

  TileFilterType::Pointer tileFilter = TileFilterType::New();

  itk::FixedArray<unsigned int, 2> layout;

  layout[0] = numberOfImgs;
  layout[1] = 0; // inputPlaneDirection.size();

  tileFilter->SetLayout( layout );
  tileFilter->SetDefaultPixelValue( 128 );
  for( unsigned int plane = 0; plane < extractingSlices.size(); plane++ )
    {
    for( unsigned int i = 0; i < numberOfImgs; i++ )
      {
      OutputRGBImageType::Pointer img = rgbSlices[i + plane * numberOfImgs];
      tileFilter->SetInput( i + plane * numberOfImgs, rgbSlices[i + plane * numberOfImgs] );
      }
    }

  /* write out 2D image */
  typedef itk::ImageFileWriter<OutputRGBImageType> RGBFileWriterType;

  RGBFileWriterType::Pointer rgbFileWriter = RGBFileWriterType::New();

  rgbFileWriter->SetInput( tileFilter->GetOutput() );
  rgbFileWriter->SetFileName( outputFilename );

  try
    {
    rgbFileWriter->Update();
    }
  catch( itk::ExceptionObject& e  )
    {
    std::cout << "ERROR:  Could not write image." << std::endl;
    std::cout << "ERROR:  " << e.what() << std::endl;
    exit(EXIT_FAILURE);
    }
  return EXIT_SUCCESS;
}
