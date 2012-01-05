#include "BRAINSCutApplyModel.h"

#include "itkLabelGeometryImageFilter.h"
#include "itkSimilarityIndexImageFilter.h"
#include "itkLabelStatisticsImageFilter.h"

#include "SimilarityIndexCLP.h"

/*
 * This is to analyse the performance of the BRAINSCut result
 * It takes in manual volume and BRAINSCut continuous volume
 * and then apply different threshold with
 * BRAINSCut's post processing method
 */

inline BinaryImagePointer ThresholdLabelImageToOneValue( BinaryImagePointer inputMaskVolume);

inline BinaryImagePointer ReadBinaryImageByFilename( std::string filename );

inline WorkingImagePointer ReadWorkingImageByFilename( std::string filename );

inline float GetVolume( BinaryImagePointer image);

void printToScreen( float manualVolume, float annVolume, float SI, float threshold);

void printHeader();

int
main(int argc, char * *argv)
{
  PARSE_ARGS;

  BRAINSCutApplyModel BRAINSCutPostProcessing;

  if( inputManualVolume == "" )
    {
    std::cout << " inputManualVolume is necessary"
              << std::endl;
    exit( EXIT_FAILURE );
    }
  /* read continuous image */
  BinaryImagePointer manualVolume = ReadBinaryImageByFilename( inputManualVolume );
  manualVolume = ThresholdLabelImageToOneValue( manualVolume );

  /* temporary file to be compared */
  BinaryImagePointer annThresholdVolume;

  /* compute manual volume */
  float floatManualVolume = GetVolume( manualVolume );

  /* set up similarity index computation */
  typedef itk::SimilarityIndexImageFilter<BinaryImageType, BinaryImageType> SimilarityIndexFilterType;
  SimilarityIndexFilterType::Pointer similarityIndexFilter = SimilarityIndexFilterType::New();

  similarityIndexFilter->SetInput1( manualVolume );

  printHeader();
  /** iterate through the threshold */
  for( float threshold = 0.0F; threshold <= 1.00F; threshold += thresholdInterval )
    {
    /* similarity index */
    annThresholdVolume = BRAINSCutPostProcessing.PostProcessingOfANNContinuousImage( ANNContinuousVolume,
                                                                                     threshold);
    similarityIndexFilter->SetInput2( annThresholdVolume );
    similarityIndexFilter->Update();

    printToScreen( floatManualVolume,
                   GetVolume( annThresholdVolume ),
                   similarityIndexFilter->GetSimilarityIndex(),
                   threshold);
    }

  return 0;
}

inline BinaryImagePointer
ThresholdLabelImageToOneValue( BinaryImagePointer inputMaskVolume)
{
  typedef itk::BinaryThresholdImageFilter<BinaryImageType, BinaryImageType> ThresholdType;
  ThresholdType::Pointer thresholder = ThresholdType::New();

  thresholder->SetInput( inputMaskVolume );
  thresholder->SetInsideValue(1);
  thresholder->SetOutsideValue(0);
  thresholder->SetLowerThreshold( 1 );
  thresholder->Update();

  BinaryImagePointer outputMask = thresholder->GetOutput();
  return outputMask;
}

inline WorkingImagePointer
ReadWorkingImageByFilename( std::string filename )
{
  typedef itk::ImageFileReader<WorkingImageType> WorkingImageReaderType;
  WorkingImageReaderType::Pointer reader = WorkingImageReaderType::New();

  reader->SetFileName( filename );
  reader->Update();

  WorkingImagePointer image = reader->GetOutput();
  return image;
}

inline BinaryImagePointer
ReadBinaryImageByFilename( std::string filename )
{
  typedef itk::ImageFileReader<BinaryImageType> BinaryImageReaderType;
  BinaryImageReaderType::Pointer reader = BinaryImageReaderType::New();

  reader->SetFileName( filename );
  reader->Update();

  BinaryImagePointer image = reader->GetOutput();
  return image;
}

inline float
GetVolume( BinaryImagePointer image)
{
  unsigned char labelValue = 1;

  typedef itk::LabelStatisticsImageFilter<BinaryImageType, BinaryImageType> MeasureFilterType;

  MeasureFilterType::Pointer manualVolumeMeasrueFilter = MeasureFilterType::New();

  manualVolumeMeasrueFilter->SetInput( image );
  manualVolumeMeasrueFilter->SetLabelInput( image );
  manualVolumeMeasrueFilter->Update();

  float count = manualVolumeMeasrueFilter->GetCount( labelValue );

  BinaryImageType::SpacingType spacing = image->GetSpacing();

  float volumeOfOneVoxel = 1.0F;
  for( unsigned int i = 0; i < DIMENSION; i++ )
    {
    volumeOfOneVoxel *= spacing[i];
    }

  return count * volumeOfOneVoxel;
}

void
printToScreen( float manualVolume,
               float annVolume,
               float SI,
               float threshold)
{
  std::cout << threshold << ", "
            << manualVolume << ", "
            << annVolume << ", "
            << SI << std::endl;
}

void
printHeader()
{
  std::cout << "threshold, manual, ann, SI" << std::endl;
}
