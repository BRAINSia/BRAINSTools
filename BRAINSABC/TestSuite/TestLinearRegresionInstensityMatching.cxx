//
// Created by Hans Johnson on 5/19/16.
//

#include "LinearRegressionIntensityMatching.h"
#include <itkImageFileReader.h>
#include <itkImage.h>


int main(int argc, char * argv[])
{
  typedef itk::Image<float, 3> ImageType;
  typedef itk::ImageFileReader<ImageType> ReaderType;
  ImageType::Pointer ref = ReaderType::New();
  ImageType::Pointer rescale =


  return EXIT_SUCCESS;
}