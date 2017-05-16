//
// Created by Hans Johnson on 5/4/17.
//
#include <algorithm>
#include <string>
#include <map>
#include <vector>
#include <fstream>
#include <array>
#include <iterator>
#include <iostream>
#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkRecursiveGaussianImageFilter.h>
#include <itkGradientImageFilter.h>
#include <itkImageRegionIteratorWithIndex.h>
#include <itkLinearInterpolateImageFunction.h>
//#include "BRAINSCutApplyModel.h"
#include "cnpy/cnpy.h"

//TODO: Write VTK lines for visualization in Slicer3D
//TODO: Fix command-line interface to use SlicerExecutionModel
//TODO: Remove hard-coded elements of this file
//TODO: Fix Gradient Magnitude mode to use "MaxGradient" of N Input feature images (or other strategy)


int main(int argc, char *argv[])
{
  if (argc > 1000000)
  {
    std::cout << "USAGE: " << argv[0]
              << "<sample_length> <sample_spacing> <label_image> <GradientMode> <first_feature_image> ... < last_feature_image>"
              << std::endl;
    //GradientMode
    //      First image, just use the first image as the gradient computation
    //      Max gradient, Magnitude.

  }
  typedef itk::Image<float, 3> FeatureImageType;
  typedef itk::Image<unsigned int, 3> LabelImageType;

  const int sample_length{std::atoi(argv[1])};
  const float sample_spacing{static_cast<float>(std::atof(argv[2]))};
  const std::string lblImageFN{argv[3]};
  const int gradient_mode{std::atoi(argv[4])};

  const int numfeatureImages{argc - 5};
  std::vector<std::string> imFeatureImageFilenames;
  for (auto i = 5; i < argc; i++)
  {
    imFeatureImageFilenames.push_back(std::string(argv[i]));
  }


  typedef itk::ImageFileReader<LabelImageType> LabelReaderType;

  LabelReaderType::Pointer labelReader = LabelReaderType::New();
  labelReader->SetFileName(lblImageFN);
  labelReader->Update();
  LabelImageType::Pointer labelImage{labelReader->GetOutput()};
  labelReader = nullptr;

  typedef itk::ImageFileReader<FeatureImageType> ReaderType;

  typedef itk::LinearInterpolateImageFunction<FeatureImageType> InterpolatorType;

  std::vector<FeatureImageType::Pointer> imVector;
  std::vector<InterpolatorType::Pointer> imInterpolatorVector;

  for (auto i = 0; i < numfeatureImages; ++i)
  {
    ReaderType::Pointer reader = ReaderType::New();
    reader->SetFileName(imFeatureImageFilenames[i]);
    reader->Update();
    FeatureImageType::Pointer image = reader->GetOutput();
    imVector.push_back(image);
    InterpolatorType::Pointer interp = InterpolatorType::New();
    interp->SetInputImage(image);
    imInterpolatorVector.push_back(interp);
  }

  typedef itk::GradientImageFilter<FeatureImageType, double, double> GradientFilterType;
  GradientFilterType::OutputImageType::Pointer gradientImage;
  if (gradient_mode == 1)
  {
    GradientFilterType::Pointer gf = GradientFilterType::New();
    gf->SetInput(imVector[0]);
    gf->Update();
    gradientImage = gf->GetOutput();
  } else
  {
    return EXIT_FAILURE;
  }

  const size_t reserve_size = labelImage->GetLargestPossibleRegion().GetNumberOfPixels() / 2;
  std::vector<float> inputVector; //continuous stream
  inputVector.reserve(reserve_size);
  std::vector<float> outputVector;
  outputVector.reserve(reserve_size);

  LabelImageType::PointType pnt;
  LabelImageType::PointType step;
  LabelImageType::PointType currPnt;

  step.Fill(1.0);
  itk::ImageRegionIteratorWithIndex<LabelImageType> it(labelImage, labelImage->GetLargestPossibleRegion());
  const float start{-(sample_spacing * (sample_length / 2))};
  const float end{-start};
  while (!it.IsAtEnd())
  {
    const float current_out_value = it.Get();

    if (current_out_value > 0.0)
    {
      LabelImageType::IndexType index = it.GetIndex();
      GradientFilterType::OutputImageType::ValueType grad = gradientImage->GetPixel(index);
      grad.Normalize();

      labelImage->TransformIndexToPhysicalPoint(index, pnt);
      for (size_t imIdx = 0; imIdx < imInterpolatorVector.size(); ++imIdx)
      {
        for (double scaleStep = start; scaleStep <= end; scaleStep += sample_spacing)
        {
          //const vnl_vector<double> temp = pnt.GetVnlVector() + scaleStep*grad.GetVnlVector();
          //currPnt.operator=(temp.data_block());
          for (size_t i = 0; i < LabelImageType::ImageDimension; ++i)
          {
            currPnt[i] = pnt[i] + scaleStep * grad[i];
          }
          float curr_value = imInterpolatorVector[imIdx]->Evaluate(currPnt);
          inputVector.push_back(curr_value);
        }
      }
      outputVector.push_back(current_out_value);
    }

    ++it;
  }

  std::cout << "InputVector:  " << inputVector.size() << std::endl;
  std::cout << "OutputVector: " << outputVector.size() << std::endl;

  unsigned int oshape[1];
  oshape[0] = outputVector.size();
  cnpy::npy_save("/Shared/johnsonhj/HDNI/20170425_CerebellumSegmentation/EditorHelperML/out.npy", &(outputVector[0]),
                 oshape, 1, "w");

  unsigned int ishape[2];
  ishape[0] = inputVector.size() / outputVector.size();
  ishape[1] = outputVector.size();
  cnpy::npy_save("/Shared/johnsonhj/HDNI/20170425_CerebellumSegmentation/EditorHelperML/in.npy", &(inputVector[0]),
                 ishape, 2, "w");

  return EXIT_SUCCESS;
}
