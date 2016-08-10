//
// Created by Johnson, Hans J on 7/29/16.
//

// \author Hans J. Johnson
// \date 2016-07-10
// Test program for evaluating matlab to ITK conversions

#include <iostream>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include "FFTWUpsample.h"

#include <itkTimeProbe.h>

#include "OpWeightedL2.h"

int main(int argc, char *argv[]) {
  if (argc != 4) {
    std::cout << "ERROR: Incorrrect number of arguments <Intensity_LR> <edgement_HR> <output>" << std::endl;
  }
  FFTWInit(""); //Use default location.

  itk::TimeProbe tp;
  tp.Start();

  const std::string highResEdgeImageFileName = argv[2];
  sitk::Image highResEdgeImage = sitk::ReadImage(highResEdgeImageFileName, sitk::sitkFloat32);

  const std::string lriFileName = argv[1];
  sitk::Image lriImage = sitk::ReadImage(lriFileName, sitk::sitkFloat32);

  sitk::Image X_lr = sitk::RescaleIntensity(lriImage, 0.0F, 1.0F);
  sitk::Image SRImage = SimpleOpWeightedL2(X_lr, highResEdgeImage);

  const std::string hriFileName = argv[3];
  sitk::Image output = sitk::Image(SRImage);
  sitk::WriteImage(output, hriFileName, sitk::sitkFloat32);

  tp.Stop();
  std::cout << tp.GetTotal() << tp.GetUnit() << std::endl;
  return EXIT_SUCCESS;
}
