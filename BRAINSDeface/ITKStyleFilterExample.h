//
// Created by Leinoff, Alexander on 5/27/16.
//

#ifndef BRAINSTOOLS_ITKSTYLEFILTEREXAMPLE_H
#define BRAINSTOOLS_ITKSTYLEFILTEREXAMPLE_H

#include "itkImage.h"
#include "itkImageToImageFilter.h"
#include <iostream>

/*
 * This is a basic example of an ITK style filter.
 * In the future, we may want to use this model for functions that
 * work within the ITK pipeline
 *
 * Example Usage:
 * #include "ITKStyleFilterExample.h"
 *
 * int main()
 * {
 *   typedef ITKStyleFilterExample <ImageType> TestFilter;
 *   TestFilter::Pointer myTest = TestFilter::New();
 *   myTest->SetInput(subject);
 *   myTest->GetOutput();
 *   myTest->Update();
 *   return 0;
 * }
 */

template<typename TInputImage >
class ITKStyleFilterExample
  : public itk::ImageToImageFilter< TInputImage, TInputImage>
{
public:
  typedef ITKStyleFilterExample Self;
  typedef itk::SmartPointer <Self> Pointer;

  itkNewMacro(Self);
  itkTypeMacro(Self, ImageToImageFilter);

protected:
  ITKStyleFilterExample(){};
  ~ITKStyleFilterExample(){};

  void GenerateData() ITK_OVERRIDE
  {
    std::cout << "Hello Filter!!!" << std::endl;
  }
};


#endif //BRAINSTOOLS_ITKSTYLEFILTEREXAMPLE_H
