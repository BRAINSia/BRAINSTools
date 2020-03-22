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
//
//
// //////////////////////////////////////////////////////////////////////////////
//
//  Registration of a dataset to an atlas using affine transformation and
//  MI image match metric
//
//  Only for 3D!
//
//  Given a list of filenames for atlas template and probabilities along with
//  the dataset, this class generate images that are in the space of the first
//  image (all data and probability images).
//
//
//
// //////////////////////////////////////////////////////////////////////////////

// prastawa@cs.unc.edu 10/2003

#ifndef __AtlasRegistrationMethod_h
#define __AtlasRegistrationMethod_h

#include "itkAffineTransform.h"
#include "itkArray.h"
#include "itkImage.h"
#include "itkObject.h"
#include "itkNaryAddImageFilter.h"

#include <vector>
#include <string>

#include "BRAINSFitHelper.h"
#include "BRAINSABCUtilities.h"
#include "itkAverageImageFilter.h"
#include "itkBinaryThresholdImageFilter.h"
#include <string>

#include "LinearRegressionIntensityMatching.h"
#include "itkIO.h"
class EmptyVectorException
{
public:
  EmptyVectorException(const char * pStr = "The list of input images was empty.  Nothing to averge.")
    : pMessage(pStr)
  {}

  const char *
  what() const
  {
    return pMessage;
  }

private:
  const char * pMessage;
};


// Take a list of coregistered images, all of the same type (T1,T2) and return the average image.
template <typename TImage>
typename TImage::Pointer
AverageImageList(const std::vector<typename TImage::Pointer> & inputImageList)
{
  if (inputImageList.empty())
  {
    // No images, something went wrong.
    throw EmptyVectorException();
  }
  if (inputImageList.size() == 1)
  {
    // Only one image, nothing to average.
    return inputImageList[0];
  }

  using BinaryThreshImageFilterType = itk::BinaryThresholdImageFilter<TImage, TImage>;
  using MultiplyFilterType = itk::MultiplyImageFilter<TImage, TImage>;
  typename BinaryThreshImageFilterType::Pointer firstBinary = BinaryThreshImageFilterType::New();
  firstBinary->SetLowerThreshold(0);
  firstBinary->SetUpperThreshold(0);
  firstBinary->SetInsideValue(0.0);
  firstBinary->SetOutsideValue(1.0);
  firstBinary->SetInput(inputImageList[0]);
  firstBinary->Update();
  typename TImage::Pointer averageMask = firstBinary->GetOutput();
  for (unsigned int i = 1; i < inputImageList.size(); ++i)
  {
    typename BinaryThreshImageFilterType::Pointer myThresholder = BinaryThreshImageFilterType::New();
    myThresholder->SetInput(inputImageList[i]);
    myThresholder->SetLowerThreshold(0); // Only valuse exactly equal to zero are to be used.
    myThresholder->SetUpperThreshold(0);
    myThresholder->SetInsideValue(0.0);
    myThresholder->SetOutsideValue(1.0);
    myThresholder->Update();
    typename MultiplyFilterType::Pointer multIF = MultiplyFilterType::New();
    multIF->SetInput1(averageMask);
    multIF->SetInput2(myThresholder->GetOutput());
    multIF->Update();
    averageMask = multIF->GetOutput();
  }

  using AvgFilterType = itk::AverageImageFilter<TImage, TImage>;
  typename AvgFilterType::Pointer filter = AvgFilterType::New();
  typename TImage::Pointer        referenceScaleImg = inputImageList[0];
  filter->SetInput(0, referenceScaleImg);
  for (unsigned int i = 1; i < inputImageList.size(); ++i)
  {
    // Modify inputImageList in place.
    typename TImage::Pointer temp = LinearRegressionIntensityMatching<TImage, TImage>(
      referenceScaleImg.GetPointer(), averageMask.GetPointer(), inputImageList[i].GetPointer());
    filter->SetInput(i, temp);
  }
  filter->Update();
  typename MultiplyFilterType::Pointer multIF = MultiplyFilterType::New();
  multIF->SetInput1(averageMask);
  multIF->SetInput2(filter->GetOutput());
  multIF->Update();

  return multIF->GetOutput();
}

/** \class AtlasRegistrationMethod
 */
template <typename TOutputPixel, typename TProbabilityPixel>
class AtlasRegistrationMethod : public itk::Object
{
public:
  /** Standard class type alias. */
  using Self = AtlasRegistrationMethod;
  using Pointer = itk::SmartPointer<Self>;
  using ConstPointer = itk::SmartPointer<const Self>;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  // Image types
  using OutputImageType = itk::Image<TOutputPixel, 3>;
  using OutputImagePointer = typename OutputImageType::Pointer;
  using OutputImageIndexType = typename OutputImageType::IndexType;
  using OutputImageOffsetType = typename OutputImageType::OffsetType;
  using OutputImagePixelType = typename OutputImageType::PixelType;
  using OutputImageSizeType = typename OutputImageType::SizeType;
  using OutputImageRegionType = typename OutputImageType::RegionType;

  using ProbabilityImageType = itk::Image<TProbabilityPixel, 3>;
  using ProbabilityImagePointer = typename ProbabilityImageType::Pointer;
  using ProbabilityImageIndexType = typename ProbabilityImageType::IndexType;
  using ProbabilityImageOffsetType = typename ProbabilityImageType::OffsetType;
  using ProbabilityImagePixelType = typename ProbabilityImageType::PixelType;
  using ProbabilityImageSizeType = typename ProbabilityImageType::SizeType;
  using ProbabilityImageRegionType = typename ProbabilityImageType::RegionType;

  using InternalImageType = itk::Image<float, 3>;
  using InternalImagePointer = typename InternalImageType::Pointer;
  using InternalImageIndexType = typename InternalImageType::IndexType;
  using InternalImageOffsetType = typename InternalImageType::OffsetType;
  using InternalImagePixelType = typename InternalImageType::PixelType;
  using InternalImageRegionType = typename InternalImageType::RegionType;
  using InternalImageSizeType = typename InternalImageType::SizeType;

  using ByteImageType = itk::Image<unsigned char, 3>;
  using ByteImagePointer = typename ByteImageType::Pointer;
  using ByteImageIndexType = typename ByteImageType::IndexType;
  using ByteImageOffsetType = typename ByteImageType::OffsetType;
  using ByteImagePixelType = typename ByteImageType::PixelType;
  using ByteImageRegionType = typename ByteImageType::RegionType;
  using ByteImageSizeType = typename ByteImageType::SizeType;

  using GenericTransformType = itk::Transform<double, 3, 3>;
  using CompositeTransformType = itk::CompositeTransform<double, 3>;
  using CompositeTransformPointer = CompositeTransformType::Pointer;
  using ProbabilityImageList = std::vector<ProbabilityImagePointer>;
  using OutputImageList = std::vector<OutputImagePointer>;

  using FlagArrayType = itk::Array<unsigned char>;

  using StringVector = std::vector<std::string>;
  using MapOfStringVectors = orderedmap<std::string, StringVector>;

  using FloatImageVector = std::vector<InternalImagePointer>;
  using MapOfFloatImageVectors = orderedmap<std::string, FloatImageVector>;

  using TransformList = std::vector<GenericTransformType::Pointer>;
  using MapOfTransformLists = orderedmap<std::string, TransformList>;

  void
  SetSuffix(std::string suffix);

  itkGetConstMacro(OutputDebugDir, std::string);
  itkSetMacro(OutputDebugDir, std::string);

  itkGetConstMacro(SaveState, std::string);
  itkSetMacro(SaveState, std::string);

  InternalImagePointer
  GetFirstAtlasOriginalImage()
  {
    return GetMapVectorFirstElement(this->m_AtlasOriginalImageList);
  }
  InternalImagePointer
  GetSecondModalityAtlasOriginalImage(const std::string & type)
  {
    MapOfFloatImageVectors::iterator test_map_location = this->m_AtlasOriginalImageList.find(type);
    if (test_map_location == this->m_AtlasOriginalImageList.end())
    {
      return nullptr;
    }
    return *(test_map_location->second.begin());
  }

  void
  SetAtlasOriginalImageList(MapOfFloatImageVectors & NewAtlasList);

  void
  SetIntraSubjectOriginalImageList(MapOfFloatImageVectors & NewIntraSubjectOriginalImageList);

  // itkSetMacro( IntraSubjectTransformFileNames, std::vector<std::string> );
  itkSetMacro(AtlasToSubjectTransformFileName, std::string);

  // INFO: KENT:  Move all code from class definition to the .hxx file outside the class definition
  void
  SetIntraSubjectTransformFileNames(MapOfStringVectors userlist)
  {
    m_IntraSubjectTransformFileNames = userlist;
    m_RegistrationUpdateNeeded = true;
  }

  void
  RegisterImages();

  GenericTransformType::Pointer
  GetAtlasToSubjectTransform()
  {
    return m_AtlasToSubjectTransform;
  }

  const MapOfTransformLists &
  GetIntraSubjectTransforms() const
  {
    return m_IntraSubjectTransforms;
  }

  // Set/Get the Debugging level for filter verboseness
  itkSetMacro(DebugLevel, unsigned int);
  itkGetMacro(DebugLevel, unsigned int);

  itkGetMacro(UseNonLinearInterpolation, bool);
  itkSetMacro(UseNonLinearInterpolation, bool);

  itkSetMacro(RestoreState, CompositeTransformPointer);
  itkGetConstMacro(RestoreState, CompositeTransformPointer);

  itkSetObjectMacro(KeySubjectImage, InternalImageType);
  itkGetModifiableObjectMacro(KeySubjectImage, InternalImageType);

  void
  SetAtlasLinearTransformChoice(const std::string & c)
  {
    m_AtlasLinearTransformChoice = c;
    m_RegistrationUpdateNeeded = true;
  }

  void
  SetImageLinearTransformChoice(const std::string & c)
  {
    m_ImageLinearTransformChoice = c;
    m_RegistrationUpdateNeeded = true;
  }

  void
  SetWarpGrid(const unsigned int gx, const unsigned int gy, const unsigned int gz)
  {
    m_WarpGrid.resize(3);
    m_WarpGrid[0] = gx;
    m_WarpGrid[1] = gy;
    m_WarpGrid[2] = gz;
    m_RegistrationUpdateNeeded = true;
  }

  void
  SetAtlasToSubjectInitialTransform(const GenericTransformType::Pointer atlasToSubjectInitialTransform)
  {
    if (this->m_AtlasToSubjectInitialTransform != atlasToSubjectInitialTransform)
    {
      this->m_AtlasToSubjectInitialTransform = atlasToSubjectInitialTransform;
      m_RegistrationUpdateNeeded = true;
    }
  }

  void
  Update();

protected:
  void
  RegisterIntraSubjectImages();
  void
  AverageIntraSubjectRegisteredImages();
  void
  RegisterAtlasToSubjectImages();

  AtlasRegistrationMethod();
  ~AtlasRegistrationMethod();

  OutputImagePointer
  CopyOutputImage(InternalImagePointer img);

  ProbabilityImagePointer
  CopyProbabilityImage(InternalImagePointer img);

private:
  std::string m_Suffix;
  std::string m_OutputDebugDir;

  //  ByteImagePointer                  m_AtlasOriginalMask;
  MapOfFloatImageVectors m_AtlasOriginalImageList;
  MapOfFloatImageVectors m_IntraSubjectOriginalImageList;
  MapOfFloatImageVectors m_RegisteredIntraSubjectImagesList;
  FloatImageVector       m_ModalityAveragedOfIntraSubjectImages;

  ByteImagePointer m_InputImageTissueRegion;
  ImageMaskPointer m_InputSpatialObjectTissueRegion;

  std::vector<unsigned int> m_WarpGrid;
  MapOfStringVectors        m_IntraSubjectTransformFileNames;
  std::string               m_AtlasToSubjectTransformFileName;

  GenericTransformType::Pointer m_AtlasToSubjectTransform;
  GenericTransformType::Pointer m_AtlasToSubjectInitialTransform;
  MapOfTransformLists           m_IntraSubjectTransforms;
  InternalImagePointer          m_KeySubjectImage; // The image to be used for intra-subject registration

  bool m_UseNonLinearInterpolation{ true };
  bool m_DoneRegistration{ false };
  bool m_RegistrationUpdateNeeded{
    true
  }; // INFO: KENT: The m_RegistrationUpdateNeeded is a hack to replicate the behavior
     // that should come from using the modified times of the itk::Object class
     //            All the Get/Set functions should use the itkSetMacro so that the
     // itk::Object->Modified times are updated correctly, then we can just use that
  //            modify status to determine when re-running is necessary.

  std::string m_AtlasLinearTransformChoice;
  std::string m_ImageLinearTransformChoice;

  std::string               m_SaveState;
  CompositeTransformPointer m_RestoreState;

  unsigned int m_DebugLevel{ 0 };
};

#ifndef MU_MANUAL_INSTANTIATION
#  include "AtlasRegistrationMethod.hxx"
#endif

#endif
