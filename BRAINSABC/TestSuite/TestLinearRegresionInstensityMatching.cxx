//
// Created by Hans Johnson on 5/19/16.
//

#include "LinearRegressionIntensityMatching.h"
#include <itkImageFileReader.h>
#include <itkImage.h>


int
main(int argc, char * argv[])
{
  using ImageType = itk::Image<float, 3>;
  using ReaderType = itk::ImageFileReader<ImageType>;
  ImageType::Pointer ref = ReaderType::New();
  ImageType::Pointer rescale =


    return EXIT_SUCCESS;
}
