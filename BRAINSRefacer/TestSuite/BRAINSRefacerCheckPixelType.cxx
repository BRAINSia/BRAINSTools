// Author: Jeffrey Obadal
/** Loads the file image information and checks if the pixel type matches that given
  *
  * Based off of Read Unknown Image type in ITKExamples
  * Latest version as of 2016/07/13: https://itk.org/ITKExamples/src/IO/ImageBase/ReadUnknownImageType/Documentation.html
  *
  *
**/

#include "BRAINSRefacerCheckPixelTypeCLP.h"
#include <itkImageIOBase.h>
#include <itkImageIOFactory.h>
#include <itkImageFileReader.h>

int main(int argc, char *argv[])
{
  PARSE_ARGS;

  std::cout << "Input image path: " << inputImageName << std::endl;
  std::cout << "Expected component type: " << pixelType << std::endl;

  itk::ImageIOBase::Pointer inputImageIO =
      itk::ImageIOFactory::CreateImageIO( inputImageName.c_str(),
                                          itk::ImageIOFactory::ReadMode );

  //TODO: add better error checking: file exists, read permissions...etc

  if( !inputImageIO )
    {
    std::cerr << "Unable to create ImageIOBase" << std::endl;
    std::cerr << "Most likely cause is the file doesn't exist" << std::endl;
    return EXIT_FAILURE;
    }

  inputImageIO->SetFileName( inputImageName );
  inputImageIO->ReadImageInformation();

  typedef itk::ImageIOBase::IOComponentType                                            IOComponentType;
  const IOComponentType realComponentType = inputImageIO->GetComponentType();
  const std::string realComponentType_string = inputImageIO->GetComponentTypeAsString( realComponentType );

  std::cout<< "Actual component type: " << realComponentType_string << std::endl;

  return realComponentType_string == pixelType ? EXIT_SUCCESS : EXIT_FAILURE;
}