#include <iostream>
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkListSample.h"
#include "itkVector.h"
#include "itkMaskImageFilter.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkCastImageFilter.h"
#include "itkFixedArray.h"
#include "itkThresholdImageFilter.h"
#include "itkImageDuplicator.h"
#include "itkMaximumRatioDecisionRule.h"
#include "itkVariableSizeMatrix.h"
#include "itkConstrainedValueDifferenceImageFilter.h"
#include "itkImageRegionConstIterator.h"
#include "itkLabelStatisticsImageFilter.h"
#include "linear.h"
#include "LogisticRegression.h"
#include "BRAINSContinuousClassCLP.h"

template <class PixelType>
int ContinuousClassification(std::string t1VolumeName, std::string T2VolumeName,
                             std::string discreteVolumeName, std::string outputVolumeName)
{
  // typedef float PixelType;
  const unsigned int Dimension = 3;

  typedef typename itk::Image<PixelType,  Dimension>    ImageType;
  typedef typename itk::Image<unsigned int, Dimension>  ShortImageType;
  typedef typename itk::ImageFileReader<ImageType>      ReaderType;
  typedef typename itk::ImageFileReader<ShortImageType> ShortReaderType;
  typedef typename itk::ImageFileWriter<ImageType>      WriterType;

  static const typename ShortImageType::PixelType grayMatterDiscreteValue = 2;
  static const typename ShortImageType::PixelType basalGrayMatterDiscreteValue = 3;
  static const typename ShortImageType::PixelType whiteMatterDiscreteValue = 1;
  static const typename ShortImageType::PixelType csfDiscreteValue = 4;
  static const typename ShortImageType::PixelType airDiscreteValue = 0;
  static const typename ShortImageType::PixelType veinousBloodDiscreteValue = 5;
  static const typename ShortImageType::PixelType allStandInDiscreteValue = 9;

  typename ReaderType::Pointer t1Reader = ReaderType::New();
  typename ReaderType::Pointer t2Reader = ReaderType::New();
  typename ShortReaderType::Pointer discreteReader = ShortReaderType::New();

  typename WriterType::Pointer outputWriter = WriterType::New();

  typename ImageType::Pointer t1Volume;
  typename ImageType::Pointer t2Volume;
  typename ShortImageType::Pointer discreteVolume;

  try
    {
    t1Reader->SetFileName( t1VolumeName );
    t1Reader->Update();
    t1Volume = t1Reader->GetOutput();

    t2Reader->SetFileName( T2VolumeName );
    t2Reader->Update();
    t2Volume = t2Reader->GetOutput();

    discreteReader->SetFileName( discreteVolumeName );
    discreteReader->Update();
    discreteVolume = discreteReader->GetOutput();
    }
  catch( itk::ExceptionObject & exe )
    {
    std::cout << exe << std::endl;
    exit(1);
    }

  // Use the labelStatistics filter to count the number of voxels for each tissue type.
  // Need this for the logistic regression problem later.

  typedef itk::LabelStatisticsImageFilter<ShortImageType, ShortImageType> LabelStatisticsFilterType;
  LabelStatisticsFilterType::Pointer labelStatisticsFilter = LabelStatisticsFilterType::New();
  labelStatisticsFilter->SetInput(discreteVolume);
  labelStatisticsFilter->SetLabelInput(discreteVolume);
  labelStatisticsFilter->Update();
  const unsigned int grayMatterSampleCount = labelStatisticsFilter->GetCount(grayMatterDiscreteValue)
    + labelStatisticsFilter->GetCount(basalGrayMatterDiscreteValue);
  const unsigned int whiteMatterSampleCount = labelStatisticsFilter->GetCount(whiteMatterDiscreteValue);
  const unsigned int csfSampleCount = labelStatisticsFilter->GetCount(csfDiscreteValue);
  const unsigned int veinousBloodSampleCount = labelStatisticsFilter->GetCount(veinousBloodDiscreteValue);

  const unsigned int            featureCount = 2; // T1 and T2
  LogisticRegression<PixelType> logisticRegressionWhiteVsCSF = LogisticRegression<PixelType>(featureCount,
                                                                                             csfSampleCount
                                                                                             + whiteMatterSampleCount);
  logisticRegressionWhiteVsCSF.SetClassOneLabel(whiteMatterDiscreteValue);
  logisticRegressionWhiteVsCSF.SetClassTwoLabel(csfDiscreteValue);
  LogisticRegression<PixelType> logisticRegressionWhiteVsGray = LogisticRegression<PixelType>(featureCount,
                                                                                              grayMatterSampleCount
                                                                                              + whiteMatterSampleCount);
  logisticRegressionWhiteVsGray.SetClassOneLabel(whiteMatterDiscreteValue);
  logisticRegressionWhiteVsGray.SetClassTwoLabel(grayMatterDiscreteValue);
  LogisticRegression<PixelType> logisticRegressionGrayVsCSF = LogisticRegression<PixelType>(featureCount,
                                                                                            grayMatterSampleCount
                                                                                            + csfSampleCount);
  logisticRegressionGrayVsCSF.SetClassOneLabel(grayMatterDiscreteValue);
  logisticRegressionGrayVsCSF.SetClassTwoLabel(csfDiscreteValue);
  LogisticRegression<PixelType> logisticRegressionVeinousBloodVsAll = LogisticRegression<PixelType>(featureCount,
                                                                                                    veinousBloodSampleCount
                                                                                                    + whiteMatterSampleCount);
  logisticRegressionVeinousBloodVsAll.SetClassOneLabel(veinousBloodDiscreteValue);
  logisticRegressionVeinousBloodVsAll.SetClassTwoLabel(allStandInDiscreteValue);

  unsigned int whiteVsGraySampleCount = 0;
  unsigned int csfVsGraySampleCount = 0;
  unsigned int whiteVsCSFSampleCount = 0;
  unsigned int veinousBloodVsAllSampleCount = 0;

  typedef itk::ImageRegionConstIterator<ImageType> ImageRegionConstIteratorType;
  ImageRegionConstIteratorType imgItr( t1Volume, t1Volume->GetRequestedRegion() );

  LogisticRegressionSample<PixelType> tempSample = LogisticRegressionSample<PixelType>(featureCount);
  std::vector<PixelType>              tempFeatures(featureCount);
  for( imgItr.GoToBegin(); !imgItr.IsAtEnd(); ++imgItr )
    {
    const typename ImageType::IndexType idx = imgItr.GetIndex();
    const typename ImageType::PixelType t1PixelValue = t1Volume->GetPixel(idx);
    const typename ImageType::PixelType t2PixelValue = t2Volume->GetPixel(idx);
    const typename ShortImageType::PixelType discretePixelValue = discreteVolume->GetPixel(idx);

    if( discretePixelValue == grayMatterDiscreteValue || discretePixelValue == basalGrayMatterDiscreteValue )
      {
      tempSample.SetLabel(grayMatterDiscreteValue);
      tempFeatures[0] = t1PixelValue;
      tempFeatures[1] = t2PixelValue;
      tempSample.SetSample(tempFeatures);
      logisticRegressionWhiteVsGray.AddLabeledSample(tempSample);

      tempSample.SetLabel(grayMatterDiscreteValue);
      tempSample.SetSample(tempFeatures);
      logisticRegressionGrayVsCSF.AddLabeledSample(tempSample);

      // tempSample.SetLabel(allStandInDiscreteValue);
      // tempSample.SetSample(tempFeatures);
      // logisticRegressionVeinousBloodVsAll.AddLabeledSample(tempSample);
      // veinousBloodVsAllSampleCount++;

      whiteVsGraySampleCount++;
      csfVsGraySampleCount++;
      }
    else if( discretePixelValue == whiteMatterDiscreteValue )
      {
      tempSample.SetLabel(whiteMatterDiscreteValue);
      tempFeatures[0] = t1PixelValue;
      tempFeatures[1] = t2PixelValue;
      tempSample.SetSample(tempFeatures);
      logisticRegressionWhiteVsGray.AddLabeledSample(tempSample);

      tempSample.SetLabel(whiteMatterDiscreteValue);
      tempSample.SetSample(tempFeatures);
      logisticRegressionWhiteVsCSF.AddLabeledSample(tempSample);

      veinousBloodVsAllSampleCount++;
      whiteVsGraySampleCount++;
      whiteVsCSFSampleCount++;
      }
    else if( discretePixelValue == csfDiscreteValue )
      {
      tempSample.SetLabel(csfDiscreteValue);
      tempFeatures[0] = t1PixelValue;
      tempFeatures[1] = t2PixelValue;
      tempSample.SetSample(tempFeatures);
      logisticRegressionGrayVsCSF.AddLabeledSample(tempSample);

      tempSample.SetLabel(csfDiscreteValue);
      tempSample.SetSample(tempFeatures);
      logisticRegressionWhiteVsCSF.AddLabeledSample(tempSample);

      // tempSample.SetLabel(allStandInDiscreteValue);
      // tempSample.SetSample(tempFeatures);
      // logisticRegressionVeinousBloodVsAll.AddLabeledSample(tempSample);
      // veinousBloodVsAllSampleCount++;

      csfVsGraySampleCount++;
      whiteVsCSFSampleCount++;
      }
    else if( discretePixelValue == veinousBloodDiscreteValue )
      {
      tempSample.SetLabel(veinousBloodDiscreteValue);
      tempFeatures[0] = t1PixelValue;
      tempFeatures[1] = t2PixelValue;
      tempSample.SetSample(tempFeatures);
      logisticRegressionVeinousBloodVsAll.AddLabeledSample(tempSample);

      veinousBloodVsAllSampleCount++;
      }
    }

  logisticRegressionWhiteVsCSF.TrainModel();
  logisticRegressionGrayVsCSF.TrainModel();
  logisticRegressionWhiteVsGray.TrainModel();
  logisticRegressionVeinousBloodVsAll.TrainModel();

  typedef typename itk::ImageDuplicator<ImageType> ImageDuplicatorType;
  typename ImageDuplicatorType::Pointer t1DuplicateImageFilter = ImageDuplicatorType::New();
  t1DuplicateImageFilter->SetInputImage(t1Reader->GetOutput() );
  t1DuplicateImageFilter->Update();
  typename ImageType::Pointer outputImage = t1DuplicateImageFilter->GetModifiableOutput();

  double predictedProbabilityEstimatesWhiteVsGray[2];
  double predictedProbabilityEstimatesWhiteVsCSF[2];
  double predictedProbabilityEstimatesGrayVsCSF[2];
  double predictedProbabilityEstimatesVeinousBloodVsAll[2];

  typename ImageType::PixelType outputAirPixelValue = 0;
  typename ImageType::PixelType outputOtherPixelValue = 9;
  typename ImageType::PixelType predictedOutputPixelValue = outputAirPixelValue;
  for( imgItr.GoToBegin(); !imgItr.IsAtEnd(); ++imgItr )
    {
    const typename ImageType::IndexType idx = imgItr.GetIndex();
    const typename ImageType::PixelType t1PixelValue = t1Reader->GetOutput()->GetPixel(idx);
    const typename ImageType::PixelType t2PixelValue = t2Reader->GetOutput()->GetPixel(idx);
    const typename ShortImageType::PixelType discretePixelValue = discreteVolume->GetPixel(idx);

    tempFeatures[0] = t1PixelValue;
    tempFeatures[1] = t2PixelValue;
    tempSample.SetSample(tempFeatures);
    logisticRegressionWhiteVsCSF.ClassifySample(tempSample);
    predictedProbabilityEstimatesWhiteVsCSF[0] = tempSample.GetLabelProbability(whiteMatterDiscreteValue);
    predictedProbabilityEstimatesWhiteVsCSF[1] = tempSample.GetLabelProbability(csfDiscreteValue);

    logisticRegressionWhiteVsGray.ClassifySample(tempSample);
    predictedProbabilityEstimatesWhiteVsGray[0] = tempSample.GetLabelProbability(whiteMatterDiscreteValue);
    predictedProbabilityEstimatesWhiteVsGray[1] = tempSample.GetLabelProbability(grayMatterDiscreteValue);

    logisticRegressionGrayVsCSF.ClassifySample(tempSample);
    predictedProbabilityEstimatesGrayVsCSF[0] = tempSample.GetLabelProbability(grayMatterDiscreteValue);
    predictedProbabilityEstimatesGrayVsCSF[1] = tempSample.GetLabelProbability(csfDiscreteValue);

    logisticRegressionVeinousBloodVsAll.ClassifySample(tempSample);
    predictedProbabilityEstimatesVeinousBloodVsAll[0] = tempSample.GetLabelProbability(veinousBloodDiscreteValue);
    predictedProbabilityEstimatesVeinousBloodVsAll[1] = tempSample.GetLabelProbability(allStandInDiscreteValue);

    if( discretePixelValue == airDiscreteValue )
      {
      predictedOutputPixelValue = outputAirPixelValue;
      }
    else if( predictedProbabilityEstimatesWhiteVsCSF[0] > predictedProbabilityEstimatesWhiteVsCSF[1] )
      {
      if( predictedProbabilityEstimatesWhiteVsGray[0] < predictedProbabilityEstimatesWhiteVsGray[1] )
        {
        if( predictedProbabilityEstimatesGrayVsCSF[0] < predictedProbabilityEstimatesGrayVsCSF[1] )
          {
          // Output voxel is other
          predictedOutputPixelValue = outputOtherPixelValue;
          }
        else
          {
          //// White is more likely, check for veinous blood?
          if( predictedProbabilityEstimatesVeinousBloodVsAll[0] > predictedProbabilityEstimatesVeinousBloodVsAll[1] )
            {
            predictedOutputPixelValue = outputOtherPixelValue;
            }
          else
            {
            // white vs gray
            predictedOutputPixelValue =
              static_cast<typename ImageType::PixelType>(130 + (120 * predictedProbabilityEstimatesWhiteVsGray[0]) );
            }
          }
        }
      else
        {
        // White is more likely, check for veinous blood?
        if( predictedProbabilityEstimatesVeinousBloodVsAll[0] < predictedProbabilityEstimatesVeinousBloodVsAll[1] )
          {
          predictedOutputPixelValue = predictedProbabilityEstimatesVeinousBloodVsAll[0] * 100;
          }
        else
          {
          // white vs gray
          predictedOutputPixelValue =
            static_cast<typename ImageType::PixelType>(130 + 120 * predictedProbabilityEstimatesWhiteVsGray[0]);
          }
        }
      }
    else
      {
      if( predictedProbabilityEstimatesGrayVsCSF[0] < predictedProbabilityEstimatesGrayVsCSF[1] )
        {
        if( predictedProbabilityEstimatesWhiteVsGray[0] > predictedProbabilityEstimatesWhiteVsGray[1] )
          {
          // Output voxel is other
          predictedOutputPixelValue = outputOtherPixelValue;
          }
        else
          {
          // CSF Vs Gray
          predictedOutputPixelValue =
            static_cast<typename ImageType::PixelType>(10 + 120 * predictedProbabilityEstimatesGrayVsCSF[0]);
          }
        }
      else
        {
        // CSF Vs Gray
        predictedOutputPixelValue =
          static_cast<typename ImageType::PixelType>(10 + 120 * predictedProbabilityEstimatesGrayVsCSF[0]);
        }
      }

    outputImage->SetPixel(idx, predictedOutputPixelValue);
    }
  std::cerr << "whiteVsGraySampleCount " << whiteVsGraySampleCount
            << " csfVsGraySampleCount " << csfVsGraySampleCount
            << " whiteVsCSFSampleCount " << whiteVsCSFSampleCount
            << " veinousBloodVsAllSampleCount " << veinousBloodVsAllSampleCount
            << std::endl;

  outputWriter->SetInput(outputImage);
  outputWriter->SetFileName(outputVolumeName);
  outputWriter->Modified();
  outputWriter->Update();

  return EXIT_SUCCESS;
}

int main(int argc, char * argv [])
{
  PARSE_ARGS;

  bool violated = false;
  if( inputT1Volume.size() == 0 )
    {
    violated = true; std::cout << "  --inputT1Volume Required! "  << std::endl;
    }
  if( inputT2Volume.size() == 0 )
    {
    violated = true; std::cout << "  --inputT2Volume Required! "  << std::endl;
    }
  if( inputDiscreteVolume.size() == 0 )
    {
    violated = true; std::cout << "  --inputDiscreteVolume Required! "  << std::endl;
    }
  if( outputVolume.size() == 0 )
    {
    violated = true; std::cout << "  --outputVolume Required! "  << std::endl;
    }
  if( violated )
    {
    exit(1);
    }

  typedef float PixelType;

  ContinuousClassification<PixelType>(inputT1Volume, inputT2Volume, inputDiscreteVolume, outputVolume);

  return EXIT_SUCCESS;
}
