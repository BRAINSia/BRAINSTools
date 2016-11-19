// Author: Jeffrey Obadal

#include "BRAINSRefacerCheckDeformationAreaLabelMapCLP.h"
#include <itkLabelImageToLabelMapFilter.h>
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkAbsoluteValueDifferenceImageFilter.h"
#include "itkStatisticsImageFilter.h"
#include "itkLabelMapMaskImageFilter.h"
#include <itkResampleImageFilter.h>
#include <itkNearestNeighborInterpolateImageFunction.h>

void outputError(itk::ExceptionObject &err)
{
  std::cerr << "Exception: " << std::endl;
  std::cerr << err << std::endl;
}

int main(int argc, char *argv[])
{
  PARSE_ARGS;

  typedef double                                                                    InputPixelType;
  const int                                                                         Dimension = 3;

  typedef itk::Image<InputPixelType, Dimension>                                     ImageType;

  typedef itk::ImageFileReader<ImageType>                                           ReaderType;

  typedef itk::LabelObject< InputPixelType, Dimension >                             LabelObjectType;
  typedef itk::LabelMap< LabelObjectType >                                          LabelMapType;
  typedef itk::LabelImageToLabelMapFilter<ImageType, LabelMapType>                  ImageToMapType;

  typedef itk::LabelMapMaskImageFilter<LabelMapType, ImageType>                     LabelMaskFilterType;

  typedef itk::AbsoluteValueDifferenceImageFilter<ImageType, ImageType, ImageType>  AbsValDiffFilterType;
  typedef itk::StatisticsImageFilter<ImageType>                                     StatisticsFilterType;

  typedef itk::NearestNeighborInterpolateImageFunction<ImageType, double>           NN_InterpolatorType;
  typedef itk::IdentityTransform<double, Dimension>                                 IdentityTransformType;
  typedef itk::ResampleImageFilter<ImageType, ImageType>                            labelResamplerType;


  // load in images
  ReaderType::Pointer originalReader = ReaderType::New();
  originalReader->SetFileName(inputOriginal);

  ReaderType::Pointer defacedReader = ReaderType::New();
  defacedReader->SetFileName(inputRefaced);

  LabelMaskFilterType::Pointer originalMaskFilter = LabelMaskFilterType::New();
  LabelMaskFilterType::Pointer defacedMaskFilter = LabelMaskFilterType::New();
  ReaderType::Pointer labelMapReader = ReaderType::New();
  ImageToMapType::Pointer imageToMapFilter = ImageToMapType::New();

  labelMapReader->SetFileName(brainLabelMap);

  //resample label map to original subject image

  NN_InterpolatorType::Pointer NN_interpolator = NN_InterpolatorType::New();
  IdentityTransformType::Pointer identityTransform = IdentityTransformType::New();
  labelResamplerType::Pointer labelResampler = labelResamplerType::New();

  labelResampler->SetInput(labelMapReader->GetOutput());
  labelResampler->SetInterpolator(NN_interpolator);
  labelResampler->SetTransform(identityTransform);
  labelResampler->SetReferenceImage(originalReader->GetOutput());
  labelResampler->UseReferenceImageOn();
  labelResampler->Update();

  // Create labelmap from label image
  imageToMapFilter->SetInput(labelResampler->GetOutput());

  //Mask the images leaving only the brain
  originalMaskFilter->SetInput(imageToMapFilter->GetOutput());
  originalMaskFilter->SetFeatureImage(originalReader->GetOutput());
  originalMaskFilter->SetLabel(0);
  if( checkNonDeformedArea )
  {
    originalMaskFilter->SetNegated(true);
  }
  originalMaskFilter->SetBackgroundValue(0);

  //Mask the images leaving only the brain
  defacedMaskFilter->SetInput(imageToMapFilter->GetOutput());
  defacedMaskFilter->SetFeatureImage(defacedReader->GetOutput());
  defacedMaskFilter->SetLabel(0);
  if( checkNonDeformedArea )
  {
    defacedMaskFilter->SetNegated(true);
  }
  defacedMaskFilter->SetBackgroundValue(0);

  AbsValDiffFilterType::Pointer absDiffFilter = AbsValDiffFilterType::New();

  absDiffFilter->SetInput1(defacedMaskFilter->GetOutput());
  absDiffFilter->SetInput2(originalMaskFilter->GetOutput());

  StatisticsFilterType::Pointer statsFilter = StatisticsFilterType::New();
  statsFilter->SetInput(absDiffFilter->GetOutput());
  try {
    statsFilter->Update();
    double absDiffSum = statsFilter->GetSum();

    std::cout << "Sum of Absolute Difference: " << absDiffSum << std::endl;

    if (checkNonDeformedArea)
    {
      if (absDiffSum == 0)
      {
        return EXIT_SUCCESS;
      }
      else
      {
        return EXIT_FAILURE;
      }
    }
    else
    {
      if( absDiffSum > 0 )
      {
        return EXIT_SUCCESS;
      }
      else
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