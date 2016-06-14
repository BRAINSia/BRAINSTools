
// Author: Jeffrey Obadal


#include "BRAINSRefacerBrainDataTestCLP.h"
#include <itkLabelImageToLabelMapFilter.h>
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkAbsoluteValueDifferenceImageFilter.h"
#include "itkStatisticsImageFilter.h"
#include "itkLabelMapMaskImageFilter.h"

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

  typedef itk::ImageFileReader<ImageType>                                           ReaderType;
  typedef itk::ImageFileWriter<ImageType>                                           WriterType;

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

  ReaderType::Pointer labelmapReader = ReaderType::New();
  labelmapReader->SetFileName(brainLabelMap);

  // Create labelmap from label image
  ImageToMapType::Pointer imageToMapFilter = ImageToMapType::New();
  imageToMapFilter->SetInput(labelmapReader->GetOutput());

  //Mask the images leaving only the brain
  LabelMaskFilterType::Pointer originalMaskFilter = LabelMaskFilterType::New();
  originalMaskFilter->SetInput(imageToMapFilter->GetOutput());
  originalMaskFilter->SetFeatureImage(originalReader->GetOutput());
  originalMaskFilter->SetLabel(0);
  originalMaskFilter->SetNegated( true );
  originalMaskFilter->SetBackgroundValue(0);

  //Mask the images leaving only the brain
  LabelMaskFilterType::Pointer defacedMaskFilter = LabelMaskFilterType::New();
  defacedMaskFilter->SetInput(imageToMapFilter->GetOutput());
  defacedMaskFilter->SetFeatureImage(defacedReader->GetOutput());
  defacedMaskFilter->SetLabel(0);
  defacedMaskFilter->SetNegated( true );
  defacedMaskFilter->SetBackgroundValue(0);

  AbsValDiffFilterType::Pointer absDiffFilter = AbsValDiffFilterType::New();
  absDiffFilter->SetInput1(defacedMaskFilter->GetOutput());
  absDiffFilter->SetInput2(originalMaskFilter->GetOutput());

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