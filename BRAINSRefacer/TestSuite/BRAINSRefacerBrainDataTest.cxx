
// Author: Jeffrey Obadal


#include "BRAINSRefacerBrainDataTestCLP.h"
#include <itkLabelImageToLabelMapFilter.h>
#include <itkMultiplyImageFilter.h>
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkAbsoluteValueDifferenceImageFilter.h"
#include "itkStatisticsImageFilter.h"
#include "itkLabelMapMaskImageFilter.h"
#include "../MaskFromLandmarksFilter.h"

void outputError(itk::ExceptionObject &err)
{
  std::cerr << "Exception: " << std::endl;
  std::cerr << err << std::endl;
}

int main(int argc, char *argv[])
{
  PARSE_ARGS;

  typedef float                                                                     InputPixelType;
  const int                                                                         Dimension = 3;

  typedef itk::Image<InputPixelType, Dimension>                                     ImageType;
  typedef itk::Image<unsigned char, Dimension >                                     MaskImageType;

  typedef itk::ImageFileReader<ImageType>                                           ReaderType;

  typedef itk::LabelObject< InputPixelType, Dimension >                             LabelObjectType;
  typedef itk::LabelMap< LabelObjectType >                                          LabelMapType;
  typedef itk::LabelImageToLabelMapFilter<ImageType, LabelMapType>                  ImageToMapType;

  typedef itk::LabelMapMaskImageFilter<LabelMapType, ImageType>                        LabelMaskFilterType;

  typedef itk::AbsoluteValueDifferenceImageFilter<ImageType, ImageType, ImageType>  AbsValDiffFilterType;
  typedef itk::StatisticsImageFilter<ImageType>                                     StatisticsFilterType;

  ReaderType::Pointer originalReader = ReaderType::New();
  originalReader->SetFileName(inputOriginal);

  ReaderType::Pointer defacedReader = ReaderType::New();
  defacedReader->SetFileName(inputRefaced);

  typedef MaskFromLandmarksFilter<ImageType, MaskImageType> MaskImageFromLandmarks;
  MaskImageFromLandmarks::Pointer landmarkReaderOriginal = MaskImageFromLandmarks::New();
  landmarkReaderOriginal->SetLandmarksFileName(brainLandmarksFile);
  landmarkReaderOriginal->SetInput(originalReader->GetOutput());

  //labelmapSwitch = true;

  MaskImageFromLandmarks::Pointer landmarkReaderDefaced = MaskImageFromLandmarks::New();
  landmarkReaderDefaced->SetLandmarksFileName(brainLandmarksFile);
  landmarkReaderDefaced->SetInput(defacedReader->GetOutput());

  //These: \/\/ are only used if using a labelmap Mask.
  LabelMaskFilterType::Pointer originalMaskFilter = LabelMaskFilterType::New();
  LabelMaskFilterType::Pointer defacedMaskFilter = LabelMaskFilterType::New();
  ReaderType::Pointer labelmapReader = ReaderType::New();
  ImageToMapType::Pointer imageToMapFilter = ImageToMapType::New();
  if( labelmapSwitch == true )
    {
    std::cout << "Using LabelMap based mask" << std::endl;
    labelmapReader->SetFileName(brainLabelMap);

    // Create labelmap from label image
    imageToMapFilter->SetInput(labelmapReader->GetOutput());

    //Mask the images leaving only the brain
    originalMaskFilter->SetInput(imageToMapFilter->GetOutput());
    originalMaskFilter->SetFeatureImage(originalReader->GetOutput());
    originalMaskFilter->SetLabel(0);
    originalMaskFilter->SetNegated(true);
    originalMaskFilter->SetBackgroundValue(0);



    //Mask the images leaving only the brain
    defacedMaskFilter->SetInput(imageToMapFilter->GetOutput());
    defacedMaskFilter->SetFeatureImage(defacedReader->GetOutput());
    defacedMaskFilter->SetLabel(0);
    defacedMaskFilter->SetNegated(true);
    defacedMaskFilter->SetBackgroundValue(0);
    }

  //multiply the images by the mask
  typedef itk::MultiplyImageFilter<MaskImageType, ImageType, ImageType > MaskMultiplyerType;
  MaskMultiplyerType::Pointer defacedMaskMultiplier = MaskMultiplyerType::New();
  defacedMaskMultiplier->SetInput1(landmarkReaderDefaced->GetOutput());
  defacedMaskMultiplier->SetInput2(defacedReader->GetOutput());

  MaskMultiplyerType::Pointer originalMaskMultiplier = MaskMultiplyerType::New();
  originalMaskMultiplier->SetInput1(landmarkReaderOriginal->GetOutput());
  originalMaskMultiplier->SetInput2(originalReader->GetOutput());

  AbsValDiffFilterType::Pointer absDiffFilter = AbsValDiffFilterType::New();

  if( labelmapSwitch == false )
    {
    std::cout << "Using landmarks based mask" << std::endl;
    absDiffFilter->SetInput1(defacedMaskMultiplier->GetOutput());
    absDiffFilter->SetInput2(originalMaskMultiplier->GetOutput());
    }
  else
    {
    absDiffFilter->SetInput1(defacedMaskFilter->GetOutput());
    absDiffFilter->SetInput2(originalMaskFilter->GetOutput());
    }


  StatisticsFilterType::Pointer statsFilter = StatisticsFilterType::New();
  statsFilter->SetInput(absDiffFilter->GetOutput());
  try
    {
    statsFilter->Update();
    double absDiffSum = statsFilter->GetSum();

    std::cout << "Sum of Absolute Difference: " << absDiffSum << std::endl;

    if( absDiffSum == 0 )
      {
      return EXIT_SUCCESS;
      }
    else
      {
      return EXIT_FAILURE;
      }

    }
  catch (itk::ExceptionObject &err)
    {
    outputError(err);
    return EXIT_FAILURE;
    }

  // should never get here
  std::cerr << "ERROR: Should never get to this point!!!" << std::endl;

  return EXIT_FAILURE;
}