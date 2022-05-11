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
#include "BRAINSMultiSTAPLECLP.h"
#include "itkIO.h"
#include "itkTransformFileReader.h"
#include "itkTransformFileWriter.h"
#include "itkIdentityTransform.h"
#include "itkResampleImageFilter.h"
#include "itkLabelImageGaussianInterpolateImageFunction.h"
#include "itkMultiLabelSTAPLEImageFilter.h"
#include "vnl/vnl_matlab_write.h"
#include <sstream>
#include <vector>

#include "BRAINSCommonLib.h"

template <typename TImage>
void
printImageStats(const TImage * image)
{
  typename TImage::RegionType region = image->GetLargestPossibleRegion();
  std::cout << "size[" << region.GetSize()[0] << "," << region.GetSize()[1] << "," << region.GetSize()[2]
            << "] Spacing[" << image->GetSpacing()[0] << "," << image->GetSpacing()[1] << "," << image->GetSpacing()[2]
            << "] Origin [" << image->GetOrigin()[0] << "," << image->GetOrigin()[1] << "," << image->GetOrigin()[2]
            << "]" << std::endl;
}

int
main(int argc, char * argv[])
{
  PARSE_ARGS;
  BRAINSRegisterAlternateIO();

  if (inputCompositeT1Volume.empty())
  {
    std::cerr << "Missing required Composite T1 Volume "
              << "use --inputCompositeT1Volume flag to specify" << std::endl;
    return 1;
  }
  if (inputLabelVolume.empty())
  {
    std::cerr << "Missing input label volumes "
              << "use --inputLabelVolume <name> to add a label volume" << std::endl;
    return 1;
  }

  if (!inputTransform.empty() && (inputLabelVolume.size() != inputTransform.size()))
  {
    std::cerr << "Transform list should have same number of" << std::endl
              << "members as as the input label volumes list" << std::endl;
    return 1;
  }

  if (outputMultiSTAPLE.empty())
  {
    std::cerr << "Missing outputMultiSTAPLE image file name" << std::endl;
    return 1;
  }

  using USImageType = itk::Image<unsigned short, 3>;

  using ImageList = std::vector<USImageType::Pointer>;
  ImageList inputLabelVolumes;
  for (const auto & it : inputLabelVolume)
  {
    USImageType::Pointer labelVolume;
    std::cout << "Reading " << it << std::endl;
    try
    {
      labelVolume = itkUtil::ReadImage<USImageType>(it);
    }
    catch (const itk::ExceptionObject & err)
    {
      std::cerr << err << std::endl;
      return 1;
    }
    inputLabelVolumes.push_back(labelVolume);
  }

  ImageList transformedLabelVolumes;

  // resample all input label images into a common space defined by
  // the input Composite volume.
  if (skipResampling)
  {
    for (const auto & l_inputLabelVolume : inputLabelVolumes)
    {
      transformedLabelVolumes.push_back(l_inputLabelVolume);
    }
  }
  else
  {
    USImageType::Pointer compositeVolume;
    try
    {
      std::cout << "Reading Composite Volume " << inputCompositeT1Volume << std::endl;
      compositeVolume = itkUtil::ReadImage<USImageType>(inputCompositeT1Volume);
    }
    catch (const itk::ExceptionObject & err)
    {
      std::cerr << err << std::endl;
      return 1;
    }
    printImageStats<USImageType>(compositeVolume);

    using TransformListType = std::vector<itk::TransformFileReader::TransformPointer>;

    TransformListType inputTransforms;

    if (!inputTransform.empty())
    {
      for (const auto & it : inputTransform)
      {
        itk::TransformFileReader::TransformPointer curTransform;

        itk::TransformFileReader::Pointer reader = itk::TransformFileReader::New();
        std::cout << "Reading " << it << std::endl;
        reader->SetFileName(it);
        try
        {
          reader->Update();
        }
        catch (const itk::ExceptionObject & err)
        {
          std::cerr << err << std::endl;
          return 1;
        }
        curTransform = reader->GetTransformList()->front();
        inputTransforms.push_back(curTransform);
      }
    }
    else
    {
      std::cout << "No transforms specified, using Identity" << std::endl;
      // fake it with identity transforms
      using IDTransformType = itk::IdentityTransform<double, 3>;

      IDTransformType::Pointer                   idXfrm = IDTransformType::New();
      itk::TransformFileReader::TransformPointer baseXfrm = idXfrm.GetPointer();
      for (std::vector<std::string>::const_iterator it = inputLabelVolume.begin(); it != inputLabelVolume.end(); ++it)
      {
        inputTransforms.push_back(baseXfrm);
      }
    }
    // set up interpolator function
    // NOTE see ANTS/Examples/make_interpolator_snip.tmp line 113 --
    // the sigma defaults to the image spacing apparently, but the
    // sigma can also be specified on the command line.
    using ucharLess = std::less<itk::NumericTraits<unsigned char>::RealType>;
    using InterpolationFunctionType = itk::LabelImageGaussianInterpolateImageFunction<USImageType, double, ucharLess>;
    InterpolationFunctionType::Pointer interpolateFunc = InterpolationFunctionType::New();
    double                             sigma[3];
    USImageType::SpacingType           spacing = compositeVolume->GetSpacing();
    for (unsigned i = 0; i < 3; ++i)
    {
      sigma[i] = spacing[i];
    }
    interpolateFunc->SetParameters(sigma, 4.0);

    std::vector<std::string>::const_iterator nameIt = inputLabelVolume.begin();
    TransformListType::const_iterator        xfrmIt = inputTransforms.begin();
    for (ImageList::const_iterator it = inputLabelVolumes.begin(); it != inputLabelVolumes.end();
         ++it, ++xfrmIt, ++nameIt)
    {
      USImageType::Pointer current = (*it);

      itk::TransformFileReader::TransformPointer curTransformBase = (*xfrmIt);

      using ResampleFilterType = itk::ResampleImageFilter<USImageType, USImageType, double>;

      std::cout << "Resampling " << (*nameIt) << std::flush;

      const auto * curTransform =
        dynamic_cast<const ResampleFilterType::TransformType *>(curTransformBase.GetPointer());
      if (curTransform == nullptr)
      {
        std::cerr << "Invalid transform " << curTransformBase << std::endl;
        exit(1);
      }
      ResampleFilterType::Pointer resampler = ResampleFilterType::New();
      try
      {
        resampler->SetInput(current);
        resampler->SetUseReferenceImage(true);
        resampler->SetReferenceImage(compositeVolume);
        resampler->SetInterpolator(interpolateFunc);
        resampler->SetTransform(curTransform);
        resampler->Update();
      }
      catch (const itk::ExceptionObject & err)
      {
        std::cerr << err << std::endl;
        return 1;
      }
      std::cout << " done." << std::endl;
      if (!resampledVolumePrefix.empty())
      {
        std::string namePart(itksys::SystemTools::GetFilenameName((*nameIt)));
        std::string resampledName = resampledVolumePrefix;
        resampledName += namePart;
        std::cerr << "Writing " << resampledName << std::flush;
        try
        {
          itkUtil::WriteImage<USImageType>(resampler->GetOutput(), resampledName);
        }
        catch (const itk::ExceptionObject & err)
        {
          std::cerr << err << std::endl;
          return 1;
        }
        std::cerr << " ... done." << std::endl;
      }
      printImageStats<USImageType>(resampler->GetOutput());
      transformedLabelVolumes.push_back(resampler->GetOutput());
    }
  }

  using STAPLEFilterType = itk::MultiLabelSTAPLEImageFilter<USImageType, USImageType>;
  STAPLEFilterType::Pointer STAPLEFilter = STAPLEFilterType::New();
  STAPLEFilter->SetNumberOfWorkUnits(1);

  if (labelForUndecidedPixels != -1)
  {
    STAPLEFilter->SetLabelForUndecidedPixels(labelForUndecidedPixels);
  }
  for (const auto & transformedLabelVolume : transformedLabelVolumes)
  {
    STAPLEFilter->PushBackInput(transformedLabelVolume);
  }

  std::cout << "Running MultiLabel Staple filter " << std::flush;
  try
  {
    STAPLEFilter->Update();
  }
  catch (const itk::ExceptionObject & err)
  {
    std::cerr << err << std::endl;
    return 1;
  }
  USImageType::Pointer output = STAPLEFilter->GetOutput();

  std::cout << " done." << std::endl;

  try
  {
    std::cout << "Writing " << outputMultiSTAPLE << std::endl;
    itkUtil::WriteImage<USImageType>(output, outputMultiSTAPLE);
  }
  catch (const itk::ExceptionObject & err)
  {
    std::cerr << err << std::endl;
    return 1;
  }

  if (!outputConfusionMatrix.empty())
  {
    std::cout << "Writing " << outputConfusionMatrix << std::endl;
    std::ofstream out(outputConfusionMatrix.c_str());
    if (!out.good())
    {
      std::cerr << "Can't write Matlab confusion matrix file " << outputConfusionMatrix << std::endl;
      return 1;
    }
    for (unsigned int i = 0; i < transformedLabelVolumes.size(); ++i)
    {
      std::stringstream name;
      name << "confusionMat" << i;
      STAPLEFilterType::ConfusionMatrixType confusionMat = STAPLEFilter->GetConfusionMatrix(i);
      vnl_matlab_write(out, confusionMat.data_array(), confusionMat.rows(), confusionMat.cols(), name.str().c_str());
    }
    out.close();
  }
  return 0;
}
