/*=========================================================================
 *
 *  Copyright SINAPSE: Scalable Informatics for Neuroscience, Processing and Software Engineering
 *            The University of Iowa
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         http://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *=========================================================================*/
#if !defined(ImageCalculatorTemplates_h__)
#  define ImageCalculatorTemplates_h__

#  define ITK_CONCEPT_NO_CHECKING
#  include "itkImageFileReader.h"
#  include "itkImageFileWriter.h"
#  include "itkImage.h"
#  include "Imgmath.h"
#  include "itkAbsImageFilter.h"
#  include "itkStatisticsImageFilter.h"
#  include "itkCastImageFilter.h"
#  include "itkOrientImageFilter.h"
#  include "itkSpatialOrientation.h"
#  include "itkMetaDataObject.h"
#  include "itkLabelStatisticsImageFilter.h"
#  include <itkSmartPointer.h>
#  include <cstdio>
#  include <cstdlib>
#  include <cstring>
#  include <sstream>
#  include <vcl_compiler.h>
#  include <iostream>
#  include "ImageCalculatorUtils.h"
#  include <metaCommand.h>
#  include <itkIdentityTransform.h>
#  include "GenericTransformImage.h"

#  define FunctorClassDeclare(name, op)                                                                                \
    template <typename PixelType>                                                                                      \
    class name                                                                                                         \
    {                                                                                                                  \
    public:                                                                                                            \
      name(const PixelType & p)                                                                                        \
        : m_Val(p){};                                                                                                  \
      name(){};                                                                                                        \
      ~name(){};                                                                                                       \
      PixelType                                                                                                        \
      operator()(const PixelType & a) const                                                                            \
      {                                                                                                                \
        return static_cast<PixelType>(a op m_Val);                                                                     \
      }                                                                                                                \
      bool                                                                                                             \
      operator==(const name & other) const                                                                             \
      {                                                                                                                \
        return this == &other;                                                                                         \
      }                                                                                                                \
      bool                                                                                                             \
      operator!=(const name & other) const                                                                             \
      {                                                                                                                \
        return !(*this == other);                                                                                      \
      }                                                                                                                \
                                                                                                                       \
    protected:                                                                                                         \
    private:                                                                                                           \
      name(const name &){};                                                                                            \
      PixelType m_Val;                                                                                                 \
    };

#  define FunctorClassDeclare2(name, op)                                                                               \
    template <typename PixelType>                                                                                      \
    class name                                                                                                         \
    {                                                                                                                  \
    public:                                                                                                            \
      name(const PixelType & p)                                                                                        \
        : m_Val(p){};                                                                                                  \
      name(){};                                                                                                        \
      ~name(){};                                                                                                       \
      PixelType                                                                                                        \
      operator()(const PixelType & a) const                                                                            \
      {                                                                                                                \
        return static_cast<PixelType>(op);                                                                             \
      }                                                                                                                \
      bool                                                                                                             \
      operator==(const name & other) const                                                                             \
      {                                                                                                                \
        return this == &other;                                                                                         \
      }                                                                                                                \
      bool                                                                                                             \
      operator!=(const name & other) const                                                                             \
      {                                                                                                                \
        return !(*this == other);                                                                                      \
      }                                                                                                                \
                                                                                                                       \
    protected:                                                                                                         \
    private:                                                                                                           \
      name(const name &){};                                                                                            \
      PixelType m_Val;                                                                                                 \
    };

namespace Functor
{
FunctorClassDeclare(mult, *);
FunctorClassDeclare(divide, /);
FunctorClassDeclare(add, +);
FunctorClassDeclare(subtract, -);
FunctorClassDeclare2(square, a * a);
FunctorClassDeclare2(binarydecimate, a > 0 ? 255 : 0);
FunctorClassDeclare2(squareroot, sqrt(static_cast<double>(a)));
} // namespace Functor

#  define FunctorProcess(op, constvalue)                                                                               \
    {                                                                                                                  \
      Functor::op<PixelType>                                                             op##functor(constvalue);      \
      typedef itk::UnaryFunctorImageFilter<ImageType, ImageType, Functor::op<PixelType>> FilterType;                   \
      typename FilterType::Pointer                                                       filter = FilterType::New();   \
      filter->SetFunctor(op##functor);                                                                                 \
      filter->SetInput(IntermediateImage);                                                                             \
      filter->Update();                                                                                                \
      IntermediateImage = filter->GetOutput();                                                                         \
    }
#  define FunctorProcess2(op)                                                                                          \
    {                                                                                                                  \
      typedef itk::UnaryFunctorImageFilter<ImageType, ImageType, Functor::op<PixelType>> FilterType;                   \
      typename FilterType::Pointer                                                       filter = FilterType::New();   \
      filter->SetInput(IntermediateImage);                                                                             \
      filter->Update();                                                                                                \
      IntermediateImage = filter->GetOutput();                                                                         \
    }

#  include "itkDiscreteGaussianImageFilter.h"
#  include "itkHistogramMatchingImageFilter.h"

static std::stringstream EffectiveInputFilters;

/*This function if called performs arithmetic operation with a constant value
 * to all the pixels in an input image.*/
template <typename ImageType>
typename ImageType::Pointer
DoGaussian(typename ImageType::Pointer input, const double sigma)
{
  using InternalImageType = itk::Image<float, ImageType::ImageDimension>;
  // Cast to float
  using ToFloatCasterType = itk::CastImageFilter<ImageType, InternalImageType>;
  typename ToFloatCasterType::Pointer toFloatCaster = ToFloatCasterType::New();
  toFloatCaster->SetInput(input);

  using FilterType = itk::DiscreteGaussianImageFilter<InternalImageType, InternalImageType>;
  /*============Filter the inputVolume using DiscreteGaussianImageFilter
   *   Include setting the x and y directions of the input images and setting order
   *     to be zero, and including normalizing Gaussian filter==================*/
  typename FilterType::Pointer filter = FilterType::New();
  filter->SetVariance(sigma);
  filter->SetMaximumError(0.01);
  filter->SetInput(toFloatCaster->GetOutput());
  filter->Update();

  using FromFloatCasterType = itk::CastImageFilter<InternalImageType, ImageType>;
  typename FromFloatCasterType::Pointer fromFloatCaster = FromFloatCasterType::New();
  fromFloatCaster->SetInput(filter->GetOutput());
  fromFloatCaster->Update();
  // Cast to data type
  return fromFloatCaster->GetOutput();
}

/*This function if called performs histogram equalization with the given number of match points*/
template <typename ImageType>
typename ImageType::Pointer
DoHisteq(typename ImageType::Pointer ref, typename ImageType::Pointer input, const int NumOfMatchPoints)
{
  using HistogramMatchingFilterType = itk::HistogramMatchingImageFilter<ImageType, ImageType>;
  typename HistogramMatchingFilterType::Pointer matchingFilter = HistogramMatchingFilterType::New();
  matchingFilter->SetSourceImage(input);
  matchingFilter->SetReferenceImage(ref);
  matchingFilter->SetNumberOfHistogramLevels(256);
  matchingFilter->SetNumberOfMatchPoints(NumOfMatchPoints);
  matchingFilter->ThresholdAtMeanIntensityOn();
  matchingFilter->Update();

  return matchingFilter->GetOutput();
}

/*This function if called performs arithmetic operation with a constant value
 * to all the pixels in an input image.*/
template <typename ImageType>
typename ImageType::Pointer
Ifilters(typename ImageType::Pointer input, MetaCommand command)
{
  using PixelType = typename ImageType::PixelType;

  typename ImageType::Pointer IntermediateImage = input;
  /*Multiplies a constant value to all the pixels of the Input image.*/
  if (!command.GetValueAsString("IMulC", "constant").empty())
  {
    const PixelType temp = static_cast<PixelType>(command.GetValueAsFloat("IMulC", "constant"));
    EffectiveInputFilters << "-ifmulc " << static_cast<double>(temp) << " ";
    FunctorProcess(mult, temp);
  }
  /*Divides a constant value from all the pixels of the Input image.*/
  if (!command.GetValueAsString("IDivC", "constant").empty())
  {
    const PixelType temp = static_cast<PixelType>(command.GetValueAsFloat("IDivC", "constant"));
    EffectiveInputFilters << "-ifdivc " << static_cast<double>(temp) << " ";
    FunctorProcess(divide, temp);
  }
  /*Adds a constant value to all the pixels of the Input image.*/
  if (!command.GetValueAsString("IAddC", "constant").empty())
  {
    const PixelType temp = static_cast<PixelType>(command.GetValueAsFloat("IAddC", "constant"));
    EffectiveInputFilters << "-ifaddc " << static_cast<double>(temp) << " ";
    FunctorProcess(add, temp);
  }

  /*Subtracts a constant value from all the pixels of the Input image.*/
  if (!command.GetValueAsString("ISubC", "constant").empty())
  {
    const PixelType temp = static_cast<PixelType>(command.GetValueAsFloat("ISubC", "constant"));
    EffectiveInputFilters << "-ifsubc " << static_cast<double>(temp) << " ";
    FunctorProcess(subtract, temp);
  }

  /*Gaussian Filters input image with value of sigma image.*/
  if (!command.GetValueAsString("IGaussianSigma", "constant").empty())
  {
    const double temp = static_cast<double>(command.GetValueAsFloat("IGaussianSigma", "constant"));
    EffectiveInputFilters << "-ifgaussiansigma " << static_cast<double>(temp) << " ";
    IntermediateImage = DoGaussian<ImageType>(IntermediateImage, temp);
  }

  /*Make Binary Output image.*/
  if (command.GetValueAsBool("Ifbin", "ifbin"))
  {
    EffectiveInputFilters << "-ifbin ";
    FunctorProcess2(binarydecimate);
  }

  /*Squares the pixels of the Input image.*/
  if (command.GetValueAsBool("ISqr", "ifsqr"))
  {
    EffectiveInputFilters << "-ifsqr ";
    FunctorProcess2(square);
  }

  /*Takes the square root of the pixels in the Input image.*/
  if (command.GetValueAsBool("ISqrt", "ifsqrt"))
  {
    EffectiveInputFilters << "-ifsqrt ";
    FunctorProcess2(squareroot);
  }
  std::cout << "--Storage type effective  input filter options:  " << EffectiveInputFilters.str() << std::endl;
  return IntermediateImage;
}

/* This function if called performs arithmetic operation with a constant value
 * to all the pixels in tee output image.*/
template <typename ImageType>
typename ImageType::Pointer
Ofilters(typename ImageType::Pointer input, MetaCommand command)
{
  using PixelType = typename ImageType::PixelType;
  std::stringstream           EffectiveOutputFilters;
  typename ImageType::Pointer IntermediateImage = input;
  /*Multiplies a constant value to all the pixels of the Output image.*/
  if (!command.GetValueAsString("OMulC", "constant").empty())
  {
    const PixelType temp = static_cast<PixelType>(command.GetValueAsFloat("OMulC", "constant"));
    EffectiveOutputFilters << "-ofmulc " << static_cast<double>(temp) << " ";
    FunctorProcess(mult, temp);
  }
  /*Divides a constant value from all the pixels of the Output image.*/
  if (!command.GetValueAsString("ODivC", "constant").empty())
  {
    const PixelType temp = static_cast<PixelType>(command.GetValueAsFloat("ODivC", "constant"));
    EffectiveOutputFilters << "-ofdivc " << static_cast<double>(temp) << " ";
    FunctorProcess(divide, temp);
  }

  /*Adds a constant value to all the pixels of the Output image.*/
  if (!command.GetValueAsString("OAddC", "constant").empty())
  {
    const PixelType temp = static_cast<PixelType>(command.GetValueAsFloat("OAddC", "constant"));
    EffectiveOutputFilters << "-ofaddc " << static_cast<double>(temp) << " ";
    FunctorProcess(add, temp);
  }

  /*Subtracts a constant value from all the pixels of the Output image.*/
  if (!command.GetValueAsString("OSubC", "constant").empty())
  {
    const PixelType temp = static_cast<PixelType>(command.GetValueAsFloat("OSubC", "constant"));
    EffectiveOutputFilters << "-ofsubc " << static_cast<double>(temp) << " ";
    FunctorProcess(subtract, temp);
  }
  /*Gaussian Filters output image with value of sigma image.*/
  if (!command.GetValueAsString("OGaussianSigma", "constant").empty())
  {
    const double temp = static_cast<double>(command.GetValueAsFloat("OGaussianSigma", "constant"));
    EffectiveOutputFilters << "-ofgaussiansigma " << static_cast<double>(temp) << " ";
    IntermediateImage = DoGaussian<ImageType>(IntermediateImage, temp);
  }

  /*Make Binary Output image.*/
  if (command.GetValueAsBool("Ofbin", "ofbin"))
  {
    EffectiveOutputFilters << "-ofbin ";
    FunctorProcess2(binarydecimate);
  }
  /*Squares the pixels of the Output image.*/
  if (command.GetValueAsBool("OSqr", "ofsqr"))
  {
    EffectiveOutputFilters << "-ofsqr ";
    FunctorProcess2(square);
  }
  /*Takes the square root of the pixels in the Output image.*/
  if (command.GetValueAsBool("OSqrt", "ofsqrt"))
  {
    EffectiveOutputFilters << "-ofsqrt ";
    FunctorProcess2(squareroot);
  }
  // return typename itk::Image< PixelType, ImageType::ImageDimension >::Pointer();
  std::cout << "--Storage type effective output filter options:  " << EffectiveOutputFilters.str() << std::endl;
  return IntermediateImage;
}

/*statfilters performs user specified statistical operations on the output image.*/
template <typename ImageType>
void
statfilters(const typename ImageType::Pointer AccImage, MetaCommand command)
{
  std::map<std::string, std::string> StatDescription{};
  std::map<std::string, float>       StatValues{};

  // The statistics image filter calclates all the statistics of AccImage
  using StatsFilterType = itk::StatisticsImageFilter<ImageType>;
  typename StatsFilterType::Pointer Statsfilter = StatsFilterType::New();
  Statsfilter->SetInput(AccImage);
  Statsfilter->Update();

  // The absolute Image filter calculates the absolute value of the pixels.
  using AbsFilterType = itk::AbsImageFilter<ImageType, ImageType>;
  typename AbsFilterType::Pointer Absfilter = AbsFilterType::New();
  Absfilter->SetInput(AccImage);
  Absfilter->Update();
  typename StatsFilterType::Pointer AbsStatsfilter = StatsFilterType::New();
  AbsStatsfilter->SetInput(Absfilter->GetOutput());
  AbsStatsfilter->Update();

  bool havestatmask = false;
  // If user gives an Input Mask Calculate the statistics of the image in the mask
  using UIntImageType = itk::Image<unsigned int, ImageType::ImageDimension>;
  using ReaderType = itk::ImageFileReader<UIntImageType>;
  typename ReaderType::Pointer reader = ReaderType::New();

  using LabelFilterType = itk::LabelStatisticsImageFilter<ImageType, UIntImageType>;
  typename LabelFilterType::Pointer MaskStatsfilter = LabelFilterType::New();
  typename LabelFilterType::Pointer MaskAbsStatsfilter = LabelFilterType::New();

  if (!command.GetValueAsString("Statmask", "File Name").empty())
  {
    reader->SetFileName(command.GetValueAsString("Statmask", "File Name").c_str());

    try
    {
      reader->Update();
    }
    catch (itk::ExceptionObject & excp)
    {
      std::cerr << "Error reading the series " << std::endl;
      std::cerr << excp << std::endl;
      throw;
    }

    havestatmask = true;

    if (command.GetValueAsString("Statmaskvalue", "constant").empty())
    {
      std::cout << "Error: If a mask image is given, a pixel value should be"
                << " entered and the Statistics in the input image will be calculated for"
                << " the pixels masked by this value.\n Skipping Statistics , Writing"
                << " output Image ." << std::endl;
      return;
    }

    MaskStatsfilter->SetInput(AccImage);
    MaskStatsfilter->SetLabelInput(reader->GetOutput());
    MaskStatsfilter->Update();

    // Absolute Values of the Masked Output Image
    typename AbsFilterType::Pointer MaskAbsfilter = AbsFilterType::New();
    MaskAbsfilter->SetInput(AccImage);
    MaskAbsfilter->Update();
    // Statistics of The Absolute Masked Output Image.
    MaskAbsStatsfilter->SetInput(MaskAbsfilter->GetOutput());
    MaskAbsStatsfilter->SetLabelInput(reader->GetOutput());
    MaskAbsStatsfilter->Update();
  }
  unsigned int MaskValue = 321231;
  if (havestatmask)
  {
    MaskValue = static_cast<unsigned int>(command.GetValueAsInt("Statmaskvalue", "constant"));
  }

  StatDescription["AVG:"] = "Average of all pixel values";
  StatDescription["MAVG:"] = "Average of all pixel values where mask > 0";

  if (command.GetValueAsBool("StatAvg", "statAVG"))
  {
    if (havestatmask)
    {
      StatValues["MAVG:"] = MaskStatsfilter->GetMean(MaskValue);
    }
    else
    {
      StatValues["AVG:"] = Statsfilter->GetMean();
    }
  }

  StatDescription["VAR:"] = "Variance of all pixel values";
  StatDescription["MVAR:"] = "Variance of all pixel values where mask > 0";

  if (command.GetValueAsBool("StatVAR", "statVAR"))
  {
    if (havestatmask)
    {
      StatValues["MVAR:"] = MaskStatsfilter->GetVariance(MaskValue);
    }
    else
    {
      StatValues["VAR:"] = Statsfilter->GetVariance();
    }
  }

  StatDescription["SUM:"] = "Sum of all pixel values";
  StatDescription["MSUM:"] = "Sum of all pixel values where mask > 0";

  if (command.GetValueAsBool("StatSUM", "statSUM"))
  {
    if (havestatmask)
    {
      StatValues["MSUM:"] = MaskStatsfilter->GetSum(MaskValue);
    }
    else
    {
      StatValues["SUM:"] = Statsfilter->GetSum();
    }
  }

  StatDescription["MIN:"] = "Minimum of all pixel values";
  StatDescription["MMIN:"] = "Minimum of all pixel values where mask > 0";

  if (command.GetValueAsBool("StatMIN", "statMIN"))
  {
    if (havestatmask)
    {
      StatValues["MMIN:"] = MaskStatsfilter->GetMinimum(MaskValue);
    }
    else
    {
      StatValues["MIN:"] = Statsfilter->GetMinimum();
    }
  }

  StatDescription["MAX:"] = "Maximum of all pixel values";
  StatDescription["MMAX:"] = "Maximum of all pixel values where mask > 0";

  if (command.GetValueAsBool("StatMAX", "statMAX"))
  {
    if (havestatmask)
    {
      StatValues["MMAX:"] = MaskStatsfilter->GetMaximum(MaskValue);
    }
    else
    {
      StatValues["MAX:"] = Statsfilter->GetMaximum();
    }
  }

  StatDescription["AMN:"] = "Minimum of the absolute value of the pixel values";
  StatDescription["MAMN:"] = "Minimum of the absolute value of the pixel values where mask > 0";

  if (command.GetValueAsBool("StatAMN", "statAMN"))
  {
    if (havestatmask)
    {
      StatValues["MAMN:"] = MaskAbsStatsfilter->GetMinimum(MaskValue);
    }
    else
    {
      StatValues["AMN:"] = AbsStatsfilter->GetMinimum();
    }
  }

  StatDescription["AMX:"] = "Maximum of the absolute value of the pixel values";
  StatDescription["MAMX:"] = "Maximum of the absolute value of the pixel values where mask > 0";

  if (command.GetValueAsBool("StatAMX", "statAMX"))
  {
    if (havestatmask)
    {
      StatValues["MAMX:"] = MaskAbsStatsfilter->GetMaximum(MaskValue);
    }
    else
    {
      StatValues["AMX:"] = AbsStatsfilter->GetMaximum();
    }
  }

  StatDescription["NPX:"] = "Number of pixels used in calculations";
  StatDescription["MNPX:"] = "Number of pixels used in calculations where mask > 0";

  if (command.GetValueAsBool("StatNPX", "statNPX"))
  {
    if (havestatmask)
    {
      typename ImageType::SizeType size;
      size = reader->GetOutput()->GetLargestPossibleRegion().GetSize();
      float NumberOfPixels = size[0] * size[1];
      if (ImageType::ImageDimension == 3)
      {
        NumberOfPixels = NumberOfPixels * size[2];
      }
      StatValues["NPX:"] = NumberOfPixels;
    }
    else
    {
      typename ImageType::SizeType size;
      size = AccImage->GetLargestPossibleRegion().GetSize();
      float NumberOfPixels = size[0] * size[1];
      if (ImageType::ImageDimension == 3)
      {
        NumberOfPixels = NumberOfPixels * size[2];
      }
      StatValues["NPX:"] = NumberOfPixels;
    }
  }
  // Show the stat values which can be calculated.
  if (command.GetValueAsBool("Statallcodes", "statallcodes"))
  {
    for (std::map<std::string, std::string>::const_iterator p = StatDescription.begin(); p != StatDescription.end();
         ++p)
    {
      std::cout << p->first << '\t' << p->second << std::endl;
    }
    std::cout << std::endl;
  }

  // Print the value map
  if ((!command.GetValueAsString("Statmask", "File Name").empty()) || command.GetValueAsBool("StatAvg", "statAVG") ||
      command.GetValueAsBool("StatVAR", "statVAR") || command.GetValueAsBool("StatSUM", "statSUM") ||
      command.GetValueAsBool("StatNPX", "statNPX"))
  {
    if (!command.GetValueAsString("OutputFilename", "filename").empty())
    {
      std::cout << "Stats for " << command.GetValueAsString("OutputFilename", "filename") << '\t';
    }
  }
  {
    for (std::map<std::string, float>::const_iterator p = StatValues.begin(); p != StatValues.end(); ++p)
    {
      std::cout << p->first << ' ' << p->second << ",  ";
    }
    std::cout << std::endl;
  }
}

/*This function is called when the user wants to write the ouput image to a file. The output image is typecasted to the
  user specified data type. */
template <typename InPixelType, typename PixelType, unsigned int ImageDims>
void
ProcessOutputStage(const typename itk::Image<InPixelType, ImageDims>::Pointer AccImage,
                   const std::string &                                        outputImageFilename,
                   MetaCommand                                                command)
{
  using InputImageType = itk::Image<InPixelType, ImageDims>;
  using OutputImageType = itk::Image<PixelType, ImageDims>;
  using CastToOutputFilterType = itk::CastImageFilter<InputImageType, OutputImageType>;
  typename CastToOutputFilterType::Pointer ToOutputTypeFilter = CastToOutputFilterType::New();

  using WriterType = itk::ImageFileWriter<OutputImageType>;
  ToOutputTypeFilter->SetInput(AccImage);
  ToOutputTypeFilter->Update();

  typename OutputImageType::Pointer OutputImage = Ofilters<OutputImageType>(ToOutputTypeFilter->GetOutput(), command);

  typename WriterType::Pointer writer = WriterType::New();
  writer->SetFileName(outputImageFilename);
  writer->SetInput(OutputImage);

  writer->Update();

  // Caluculate Statistics of the Image.
  statfilters<OutputImageType>(OutputImage, command);
}

class string_tokenizer : public std::vector<std::string>
{
public:
  explicit string_tokenizer(const std::string & s, const char * const sep = " ") { this->init(s, sep); }
  string_tokenizer & operator=(const string_tokenizer &)  = delete;
  string_tokenizer(const string_tokenizer &) = delete;

  explicit string_tokenizer(const char * const s, const char * const sep = " ") { this->init(std::string(s), sep); }
private:
  void
  init(const std::string & input, const char * const sep = " ")
  {
    std::string::size_type start = 0;

    std::string::size_type _end = 0;
    int                    i = 0;

    while ((start = input.find_first_not_of(sep, _end)) != std::string::npos)
    {
      this->resize(i + 1);
      _end = input.find_first_of(sep, start);
      (*this)[i] = input.substr(start, _end - start);
      ++i;
    }
  }
};

template <typename ImageType>
typename ImageType::Pointer
DoResampleStep(int                         interpCode,
               typename ImageType::Pointer ref_space,
               bool                        isbinary,
               typename ImageType::Pointer imageToResample)
{
  static typename itk::IdentityTransform<double, ImageType::ImageDimension>::Pointer id_tfm =
    itk::IdentityTransform<double, ImageType::ImageDimension>::New();
  if (ref_space.IsNull())
  {
    return imageToResample;
  }

  /*If the accumulator buffer is not empty, then every subsequent image is resampled to the current
  accumulator buffer.*/
  std::string interpolationMode{ "UNKOWNN" };
  switch (interpCode)
  {
    case 0:
      interpolationMode = "NearestNeighbor";
      break;
    case 1:
      interpolationMode = "Linear";
      break;
    case 2:
      interpolationMode = "WindowedSinc";
      break;
    default:
      itkGenericExceptionMacro("Unkown interpolation mode specified");
  }

  ImageType * const inimgptr = imageToResample.GetPointer();
  ImageType * const refimgptr = ref_space.GetPointer();
  typename itk::Transform<double, ImageType::ImageDimension, ImageType::ImageDimension>::ConstPointer defaulttfm =
    dynamic_cast<typename itk::Transform<double, ImageType::ImageDimension, ImageType::ImageDimension> *>(
      id_tfm.GetPointer());
  imageToResample =
    SimpleGenericTransformImage<ImageType>(inimgptr, refimgptr, defaulttfm, 0, interpolationMode, isbinary);
  return imageToResample;
}


/*This function reads in the input images and writes the output image ,
 * delegating the computations to other functions*/
template <typename ImageType>
void
ImageCalculatorReadWrite(MetaCommand & command)
{
  // Replace backslash blank with a unique string
  std::string tempFilenames = command.GetValueAsString("in");


  const bool isbinary = command.GetValueAsBool("Ifbin", "ifbin");

  const bool interp_mode_set = command.GetOptionWasSet("IResample");
  const int  interpCode = (interp_mode_set) ? static_cast<int>(command.GetValueAsInt("IResample", "constant")) : -1;
  if (interp_mode_set)
  {
    EffectiveInputFilters << "-ifresample " << static_cast<int>(interpCode) << " ";
  }

  const bool        reference_resample_set = command.GetOptionWasSet("IResampleRef");
  const std::string ref_space_filename =
    (reference_resample_set) ? command.GetValueAsString("IResampleRef", "constant") : "";
  if (reference_resample_set)
  {
    EffectiveInputFilters << "-ifresampleref " << ref_space_filename << " ";
  }
  typename ImageType::Pointer ref_space =
    (reference_resample_set) ? itkUtil::ReadImage<ImageType>(ref_space_filename) : nullptr;

  ReplaceSubWithSub(tempFilenames, "\\ ", "BACKSLASH_BLANK");
  // Now split into separate filenames
  string_tokenizer InputList(tempFilenames, " ");
  // Finally, return the blanks to the filenames
  for (size_t i = 0; i < InputList.size(); ++i)
  {
    ReplaceSubWithSub(InputList[i], "BACKSLASH_BLANK", " ");
  }

  using ReaderType = itk::ImageFileReader<ImageType>;
  using PixelType = typename ImageType::PixelType;
  // Read the first Image
  typename ReaderType::Pointer reader= ReaderType::New();
  reader->SetFileName(InputList.at(0).c_str());
  std::cout << "Reading 1st Image..." << InputList.at(0).c_str() << std::endl;
  try
  {
    reader->Update();
  }
  catch (itk::ExceptionObject & excp)
  {
    std::cerr << "Error reading the series " << excp << std::endl;
    throw;
  }
  // Create an Accumulator Image.
  typename ImageType::Pointer AccImage = Ifilters<ImageType>(reader->GetOutput(), command);
  AccImage = DoResampleStep<ImageType>(interpCode, ref_space, isbinary, AccImage);

  /*For variance image first step is to square the input image.*/
  typename ImageType::Pointer SqrImageSum;
  if (command.GetValueAsBool("Var", "var"))
  {
    SqrImageSum = Imul<ImageType>(reader->GetOutput(), reader->GetOutput());
  }
  /* Accumulator contains the first image initially and is updated by the next image at each count */
  for (unsigned int currimage = 1; currimage < InputList.size(); ++currimage)
  {
    std::cout << "Reading image.... " << InputList.at(currimage).c_str() << std::endl;

    typename ReaderType::Pointer reader2= ReaderType::New();
    reader2->SetFileName(InputList.at(currimage).c_str());
    try
    {
      reader2->Update();
    }
    catch (itk::ExceptionObject & excp)
    {
      std::cerr << "Error reading the series " << excp << std::endl;
      throw;
    }
    typename ImageType::Pointer SubSequentImage = reader2->GetOutput();
    reader2 = nullptr;

    /*If the accumulator buffer is not empty, then every subsequent image is histogram equalized to the current
      accumulator buffer.*/
    if (!command.GetValueAsString("IHisteq", "constant").empty() && AccImage.IsNotNull())
    {
      const int NumOfMatchPoints = static_cast<int>(command.GetValueAsInt("IHisteq", "constant"));
      EffectiveInputFilters << "-ifhisteq " << static_cast<int>(NumOfMatchPoints) << " ";
      SubSequentImage = DoHisteq<ImageType>(AccImage, SubSequentImage, NumOfMatchPoints);
    }
    SubSequentImage = DoResampleStep<ImageType>(interpCode, ref_space, isbinary, SubSequentImage);

    // Check whether the image dimensions and the spacing are the same.
    if ((AccImage->GetLargestPossibleRegion().GetSize() != SubSequentImage->GetLargestPossibleRegion().GetSize()))
    {
      itkGenericExceptionMacro(<< "Error:: The size of the images don't match.");
    }

    typename ImageType::Pointer image = Ifilters<ImageType>(SubSequentImage, command);

    vnl_vector_fixed<double, 3> spacingDifference;
    spacingDifference[0] = AccImage->GetSpacing()[0] - image->GetSpacing()[0];
    spacingDifference[1] = AccImage->GetSpacing()[1] - image->GetSpacing()[1];
    spacingDifference[2] = AccImage->GetSpacing()[2] - image->GetSpacing()[2];

    if (spacingDifference.two_norm() > 0.0001) // HACK:  Should be a percentage of the actual spacing size.
    {
      itkGenericExceptionMacro(<< "ERROR:: The pixel spacing of the images are not close enough.");
    }
    else if (AccImage->GetSpacing() != image->GetSpacing())
    {
      std::cout << "WARNING: ::The pixel spacing of the images don't match exactly. \n";
    }
    if (AccImage->GetDirection() != image->GetDirection())
    {
      itkGenericExceptionMacro(<< "Error:: The orientation of the images are different.");
    }

    // Do the math for the Accumulator image and the image read in for each iteration.
    /*Call the multiplication function*/
    if (command.GetValueAsBool("Mul", "mul"))
    {
      AccImage = Imul<ImageType>(AccImage, image);
    }

    /*Call the addition function*/
    if (command.GetValueAsBool("Add", "add"))
    {
      AccImage = Iadd<ImageType>(AccImage, image);
    }

    /*Call the subtraction function*/
    if (command.GetValueAsBool("Sub", "sub"))
    {
      AccImage = Isub<ImageType>(AccImage, image);
    }

    /*Call the divisionion function*/
    if (command.GetValueAsBool("Div", "div"))
    {
      AccImage = Idiv<ImageType>(AccImage, image);
    }

    /*For Average we add the images first*/
    if (command.GetValueAsBool("Avg", "avg"))
    {
      AccImage = Iadd<ImageType>(AccImage, image);
    }

    /*For variance we add the square image to the image in the current
      iteration and store the sum of the square image.*/
    if (command.GetValueAsBool("Var", "var"))
    {
      AccImage = Iadd<ImageType>(AccImage, image);
      typename ImageType::Pointer multimage = Imul<ImageType>(image, image);
      SqrImageSum = Iadd<ImageType>(SqrImageSum, multimage);
    }
  }

  const int NumImages = InputList.size();
  // To get the average image we divide the accumulator image with the total number of images.
  if (command.GetValueAsBool("Avg", "avg"))
  {
    AccImage = Iavg<ImageType>(AccImage, NumImages);
  }

  // Image variance is calculated.
  if (command.GetValueAsBool("Var", "var"))
  {
    typename ImageType::Pointer NumSqrImageSum =
      ImageMultiplyConstant<ImageType>(SqrImageSum, static_cast<PixelType>(NumImages));
    AccImage = Imul<ImageType>(AccImage, AccImage);
    AccImage = Isub<ImageType>(NumSqrImageSum, AccImage);
    AccImage = ImageDivideConstant<ImageType>(AccImage, static_cast<PixelType>(NumImages * NumImages - NumImages));
  }

  const std::string OutType(command.GetValueAsString("OutputPixelType", "PixelType"));

  // The resultant Image is written.
  if (!command.GetValueAsString("OutputFilename", "filename").empty())
  {
    const std::string outputFilename(command.GetValueAsString("OutputFilename", "filename"));
    std::cout << "Before write..." << outputFilename << std::endl;
    // Type cast image according to the output type specified. Default is float.
    if (!command.GetValueAsString("OutputPixelType", "PixelType").empty())
    {
      // process the string for the data type
      if (CompareNoCase(OutType, std::string("UCHAR")) == 0)
      {
        ProcessOutputStage<PixelType, unsigned char, ImageType::ImageDimension>(AccImage, outputFilename, command);
      }
      else if (CompareNoCase(OutType, std::string("SHORT")) == 0)
      {
        ProcessOutputStage<PixelType, short, ImageType::ImageDimension>(AccImage, outputFilename, command);
      }
      else if (CompareNoCase(OutType, std::string("USHORT")) == 0)
      {
        ProcessOutputStage<PixelType, unsigned short, ImageType::ImageDimension>(AccImage, outputFilename, command);
      }
      else if (CompareNoCase(OutType, std::string("INT")) == 0)
      {
        ProcessOutputStage<PixelType, int, ImageType::ImageDimension>(AccImage, outputFilename, command);
      }
      else if (CompareNoCase(OutType, std::string("UINT")) == 0)
      {
        ProcessOutputStage<PixelType, unsigned int, ImageType::ImageDimension>(AccImage, outputFilename, command);
      }
      else if (CompareNoCase(OutType, std::string("FLOAT")) == 0)
      {
        ProcessOutputStage<PixelType, float, ImageType::ImageDimension>(AccImage, outputFilename, command);
      }
      else if (CompareNoCase(OutType, std::string("DOUBLE")) == 0)
      {
        ProcessOutputStage<PixelType, double, ImageType::ImageDimension>(AccImage, outputFilename, command);
      }
      else
      {
        std::cout << "Error. Invalid data type for -outtype!  Use one of these:" << std::endl;
        PrintDataTypeStrings();
        throw;
      }
    }
    else
    {
      // Default is the Input Pixel Type.
      ProcessOutputStage<PixelType, PixelType, ImageType::ImageDimension>(AccImage, outputFilename, command);
    }
  }
}

/*This function calls the ImageCalculatorReadWrite function based on the data type specified by the user.*/
template <unsigned int DIMS>
void
ImageCalculatorProcessND(const std::string & InType, MetaCommand & command)
{
  if (CompareNoCase(InType, std::string("UCHAR")) == 0)
  {
    ImageCalculatorReadWrite<itk::Image<unsigned char, DIMS>>(command);
  }
  else if (CompareNoCase(InType, std::string("SHORT")) == 0)
  {
    ImageCalculatorReadWrite<itk::Image<short, DIMS>>(command);
  }
  else if (CompareNoCase(InType, std::string("USHORT")) == 0)
  {
    ImageCalculatorReadWrite<itk::Image<unsigned short, DIMS>>(command);
  }
  else if (CompareNoCase(InType, std::string("INT")) == 0)
  {
    ImageCalculatorReadWrite<itk::Image<int, DIMS>>(command);
  }
  else if (CompareNoCase(InType, std::string("UINT")) == 0)
  {
    ImageCalculatorReadWrite<itk::Image<unsigned int, DIMS>>(command);
  }
  else if (CompareNoCase(InType, std::string("FLOAT")) == 0)
  {
    ImageCalculatorReadWrite<itk::Image<float, DIMS>>(command);
  }
  else if (CompareNoCase(InType, std::string("DOUBLE")) == 0)
  {
    ImageCalculatorReadWrite<itk::Image<double, DIMS>>(command);
  }
}

#endif // ImageCalculatorTemplates_h__
