// Author: Jeffrey Obadal


#include "BRAINSRefacerCheckDeformationAreaLandmarksCLP.h"
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

int main(int argc, char *argv[]) {
  PARSE_ARGS;

  const bool checkDeformedArea = !checkNonDeformedArea;

  typedef float                                                                           InputPixelType;
  const int                                                                               Dimension = 3;

  typedef itk::Image<InputPixelType, Dimension>                                           ImageType;
  typedef itk::Image<unsigned char, Dimension>                                            MaskImageType;

  typedef itk::ImageFileReader<ImageType>                                                 ReaderType;

  typedef itk::AbsoluteValueDifferenceImageFilter<ImageType, ImageType, ImageType>        AbsValDiffFilterType;
  typedef itk::StatisticsImageFilter<ImageType> StatisticsFilterType;

  ReaderType::Pointer originalReader = ReaderType::New();
  originalReader->SetFileName(inputOriginal);

  ReaderType::Pointer defacedReader = ReaderType::New();
  defacedReader->SetFileName(inputRefaced);

  typedef MaskFromLandmarksFilter<ImageType, MaskImageType> MaskImageFromLandmarksFilterType;
  typedef itk::MultiplyImageFilter<MaskImageType, ImageType, ImageType> MaskMultiplierType;

  // setup original image masked for deformation area
  MaskImageFromLandmarksFilterType::Pointer areaMaskFilterOriginal = MaskImageFromLandmarksFilterType::New();
  areaMaskFilterOriginal->SetLandmarksFileName(brainLandmarksFile);
  areaMaskFilterOriginal->SetInput(originalReader->GetOutput());

  MaskMultiplierType::Pointer areaOriginalMaskMultiplier = MaskMultiplierType::New();
  areaOriginalMaskMultiplier->SetInput1(areaMaskFilterOriginal->GetOutput());
  areaOriginalMaskMultiplier->SetInput2(originalReader->GetOutput());

  // setup refaced image masked for deformation area
  MaskImageFromLandmarksFilterType::Pointer areaMaskFilterRefaced = MaskImageFromLandmarksFilterType::New();
  areaMaskFilterRefaced->SetLandmarksFileName(brainLandmarksFile);
  areaMaskFilterRefaced->SetInput(defacedReader->GetOutput());

  //set which area to check
  if (checkDeformedArea)
  {
    areaMaskFilterOriginal->SetReverseMask(true);
    areaMaskFilterRefaced->SetReverseMask(true);
  }

  MaskMultiplierType::Pointer areaRefacedMaskMultiplier = MaskMultiplierType::New();
  areaRefacedMaskMultiplier->SetInput1(areaMaskFilterRefaced->GetOutput());
  areaRefacedMaskMultiplier->SetInput2(defacedReader->GetOutput());

  // Get the absolute value difference of the two deformed area masked images
  AbsValDiffFilterType::Pointer absDiffAreaFilter = AbsValDiffFilterType::New();
  absDiffAreaFilter->SetInput1(areaOriginalMaskMultiplier->GetOutput());
  absDiffAreaFilter->SetInput2(areaRefacedMaskMultiplier->GetOutput());

  // Setup the statistics filter
  StatisticsFilterType::Pointer areaStatsFilter = StatisticsFilterType::New();
  areaStatsFilter->SetInput(absDiffAreaFilter->GetOutput());

  try
  {
    areaStatsFilter->Update();
    double absDiffSum = areaStatsFilter->GetSum();

    std::cout << "Sum of Absolute Difference: " << absDiffSum << std::endl;

    if( checkDeformedArea )
    {
      //
      return absDiffSum > 0 ? EXIT_SUCCESS : EXIT_FAILURE;
    }

    else // check non deformed area
    {
      // the brain is in this area, and thus there should be no difference
      return absDiffSum == 0 ? EXIT_SUCCESS : EXIT_FAILURE;
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