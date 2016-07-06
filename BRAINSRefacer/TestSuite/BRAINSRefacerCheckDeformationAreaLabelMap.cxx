// Author: Jeffrey Obadal

#include "BRAINSRefacerCheckDeformationAreaLabelMapCLP.h"
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

  typedef itk::LabelObject< InputPixelType, Dimension >                             LabelObjectType;
  typedef itk::LabelMap< LabelObjectType >                                          LabelMapType;
  typedef itk::LabelImageToLabelMapFilter<ImageType, LabelMapType>                  ImageToMapType;

  typedef itk::LabelMapMaskImageFilter<LabelMapType, ImageType>                     LabelMaskFilterType;

  typedef itk::AbsoluteValueDifferenceImageFilter<ImageType, ImageType, ImageType>  AbsValDiffFilterType;
  typedef itk::StatisticsImageFilter<ImageType>                                     StatisticsFilterType;

  ReaderType::Pointer originalReader = ReaderType::New();
  originalReader->SetFileName(inputOriginal);

  ReaderType::Pointer defacedReader = ReaderType::New();
  defacedReader->SetFileName(inputRefaced);

  LabelMaskFilterType::Pointer originalMaskFilter = LabelMaskFilterType::New();
  LabelMaskFilterType::Pointer defacedMaskFilter = LabelMaskFilterType::New();

  ReaderType::Pointer labelMapReader = ReaderType::New();
  ImageToMapType::Pointer imageToMapFilter = ImageToMapType::New();

  labelMapReader->SetFileName(brainLabelMap);

  // Create label map from label image
  imageToMapFilter->SetInput(labelMapReader->GetOutput());

  //Mask the images
  //Note that this actually masks the non-brain data
  //which is the deformation area

  originalMaskFilter->SetInput(imageToMapFilter->GetOutput());
  originalMaskFilter->SetFeatureImage(originalReader->GetOutput());
  originalMaskFilter->SetLabel(0);
  originalMaskFilter->SetBackgroundValue(0);

  defacedMaskFilter->SetInput(imageToMapFilter->GetOutput());
  defacedMaskFilter->SetFeatureImage(defacedReader->GetOutput());
  defacedMaskFilter->SetLabel(0);
  defacedMaskFilter->SetBackgroundValue(0);

  // reverse the mask for checking the brain data
  if( checkNonDeformedArea )
    {
    originalMaskFilter->SetNegated(true);
    defacedMaskFilter->SetNegated(true);
    }


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


    if( checkNonDeformedArea )
      {
      //NonDeformedArea includes the brain, so there should be zero differences
      return (absDiffSum == 0 ) ? EXIT_SUCCESS : EXIT_FAILURE;
      }
    else //check deformed area
      {
      // We want some deformation here, so if it's doing something it should be > 0
      return ( absDiffSum > 0 ) ? EXIT_SUCCESS : EXIT_FAILURE;
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