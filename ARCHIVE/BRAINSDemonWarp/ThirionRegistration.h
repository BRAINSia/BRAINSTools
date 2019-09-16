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
#ifndef __ThirionRegistration_h
#define __ThirionRegistration_h

#include <string>
#include "DemonsPreprocessor.h"
#include "DemonsRegistrator.h"
#include "ValidationInputParser.h"
#include "ApplicationBase.h"
#include "itkCheckerBoardImageFilter.h"
namespace itk
{
/*This file defines Thirion registration class which initializes the input
 * parser, preprocessor and the registrator. */

template <typename TImage, typename TRealImage, typename TOutputImage>
class ThirionRegistration
  : public ApplicationBase<ValidationInputParser<TImage>,
                           DemonsPreprocessor<TImage, TRealImage>,
                           DemonsRegistrator<TRealImage, TOutputImage, typename TRealImage::PixelType>>
{
public:
  /** Standard class type alias. */
  using Self = ThirionRegistration;
  using Superclass = ApplicationBase<ValidationInputParser<TImage>,
                                     DemonsPreprocessor<TImage, TRealImage>,
                                     DemonsRegistrator<TRealImage, TRealImage, typename TRealImage::PixelType>>;
  using Pointer = SmartPointer<Self>;
  using ConstPointer = SmartPointer<const Self>;

  /** Displacement field pixel type. */
  using FieldValueType = float;
  using FieldPixelType = Vector<FieldValueType, 3>;
  using TDisplacementField = Image<FieldPixelType, 3>;

  /** Run-time type information (and related methods). */
  itkTypeMacro(ThirionRegistration, ApplicationBase);

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Image types. */
  using ImageType = TImage;
  using RealImageType = TRealImage;

  /** Image dimension. */
  static constexpr unsigned int ImageDimension = TImage::ImageDimension;

  /** Type to hold the number of checker boxes per dimension */
  using PatternArrayType = FixedArray<unsigned int, TImage::ImageDimension>;

  using PixelType = typename ImageType::PixelType;
  using IndexType = typename ImageType::IndexType;
  using SizeType = typename ImageType::SizeType;

  /** ShrinkFactors type. */
  using ShrinkFactorsType = FixedArray<unsigned int, TImage::ImageDimension>;

  /** IterationArray type. */
  using IterationsArrayType = Array<unsigned int>;

  /** Set the atlas patient ID. */
  itkSetStringMacro(TheMovingImageFilename);
  itkGetStringMacro(TheMovingImageFilename);

  /** Set the subject patient ID. */
  itkSetStringMacro(TheFixedImageFilename);
  itkGetStringMacro(TheFixedImageFilename);

  /** Set the initial Displacement Field one of 3 ways. */
  itkSetStringMacro(InitialDisplacementFieldFilename);
  itkGetStringMacro(InitialDisplacementFieldFilename);

  itkSetStringMacro(InitialCoefficientFilename);
  itkGetStringMacro(InitialCoefficientFilename);

  itkSetStringMacro(InitialTransformFilename);
  itkGetStringMacro(InitialTransformFilename);

  /** Set Displacementname */
  itkSetStringMacro(DisplacementBaseName);
  itkGetStringMacro(DisplacementBaseName);
  /** Set WarpedImageName */
  itkSetStringMacro(WarpedImageName);
  itkGetStringMacro(WarpedImageName);

  /** Set input parameter file */
  // itkSetStringMacro (ParameterFilename);

  /** Set output transformation filename. */
  itkSetStringMacro(OutputFilename);

  /** Set checker board Image filename */
  itkSetStringMacro(CheckerBoardFilename);
  itkGetStringMacro(CheckerBoardFilename);

  /** Set Displacement field output filename */
  itkSetStringMacro(DisplacementFieldOutputName);
  itkGetStringMacro(DisplacementFieldOutputName);

  /** Set Checker pattern */
  itkSetMacro(CheckerBoardPattern, PatternArrayType);
  itkGetConstReferenceMacro(CheckerBoardPattern, PatternArrayType);

  /** Set append output file boolean. */
  itkSetMacro(AppendOutputFile, bool);
  itkGetMacro(AppendOutputFile, bool);
  itkBooleanMacro(AppendOutputFile);

  /* BOBF macros
   * Set Target Mask filename */
  itkSetStringMacro(BOBFTargetMask);
  itkGetStringMacro(BOBFTargetMask);

  /** Set Template Mask filename */
  itkSetStringMacro(BOBFTemplateMask);
  itkGetStringMacro(BOBFTemplateMask);

  /** Force Centered Image. */
  itkSetMacro(ForceCoronalZeroOrigin, bool);
  itkGetConstMacro(ForceCoronalZeroOrigin, bool);

  /** Output Normalized Image. */
  itkSetStringMacro(OutNormalized);
  itkGetStringMacro(OutNormalized);

  /** Set/Get the lower threshold. The default is 0. */
  itkSetMacro(Lower, PixelType);
  itkGetMacro(Lower, PixelType);

  /** Set/Get the upper threshold. The default is 70 */
  itkSetMacro(Upper, PixelType);
  itkGetMacro(Upper, PixelType);

  /** Set/Get value to replace thresholded pixels. Pixels that lie *
   *  within Lower and Upper (inclusive) will be replaced with this
   *  value. The default is 1. */
  itkSetMacro(DefaultPixelValue, PixelType);
  itkGetMacro(DefaultPixelValue, PixelType);

  /** Set the radius of the neighborhood used for a mask. */
  itkSetMacro(Radius, SizeType);
  /** Get the radius of the neighborhood used to compute the median */
  itkGetConstReferenceMacro(Radius, SizeType);

  /** Set the Seed of the neighborhood used for a mask. */
  itkSetMacro(Seed, IndexType);
  /** Get the radius of the neighborhood used to compute the median */
  itkGetConstReferenceMacro(Seed, IndexType);

  itkSetMacro(MedianFilterSize, SizeType);
  itkGetMacro(MedianFilterSize, SizeType);

  /** Set the initial deformation field to prime registration */
  //    itkSetObjectMacro(InitialDisplacementField,TDisplacementField);
  /** Set the Input Landmark Filename*/
  itkSetStringMacro(FixedLandmarkFilename);
  itkGetStringMacro(FixedLandmarkFilename);
  itkSetStringMacro(MovingLandmarkFilename);
  itkGetStringMacro(MovingLandmarkFilename);

  /** Set histogram matching*/
  itkSetMacro(UseHistogramMatching, bool);
  itkGetConstMacro(UseHistogramMatching, bool);

  /** Get the number of histogram bins. */
  itkGetConstMacro(NumberOfHistogramLevels, unsigned long);
  itkSetMacro(NumberOfHistogramLevels, unsigned long);

  /** Get the number of match points. */
  itkGetConstMacro(NumberOfMatchPoints, unsigned long);
  itkSetMacro(NumberOfMatchPoints, unsigned long);

  /** Get the number of levels. */
  itkGetMacro(NumberOfLevels, unsigned short);
  itkSetMacro(NumberOfLevels, unsigned short);

  /** Get the atlas image starting shrink factors. */
  itkGetConstReferenceMacro(TheMovingImageShrinkFactors, ShrinkFactorsType);
  void
  SetTheMovingImageShrinkFactors(const ShrinkFactorsType & shrinkfactors)
  {
    this->m_TheMovingImageShrinkFactors = shrinkfactors;
  }

  /** Get the subject image starting shrink factors. */
  itkGetConstReferenceMacro(TheFixedImageShrinkFactors, ShrinkFactorsType);
  void
  SetTheFixedImageShrinkFactors(const ShrinkFactorsType & shrinkfactors)
  {
    this->m_TheFixedImageShrinkFactors = shrinkfactors;
  }

  /** Get the number of iterations at each level. */
  itkGetConstReferenceMacro(NumberOfIterations, IterationsArrayType);
  void
  SetNumberOfIterations(const IterationsArrayType & iterations)
  {
    m_NumberOfIterations = iterations;
  }

  using BaseRegistrationFilterType =
    itk::PDEDeformableRegistrationFilter<RealImageType, RealImageType, TDisplacementField>;
  void
  SetRegistrationFilter(typename BaseRegistrationFilterType::Pointer filter)
  {
    this->m_Registrator->SetRegistrationFilter(filter);
  }

protected:
  ThirionRegistration();
  virtual ~ThirionRegistration() {}

  /** Initialize the input parser. */
  virtual void
  InitializeParser();

  /** Initialize the preprocessor */
  virtual void
  InitializePreprocessor();

  /** Initialize the registrator  */
  virtual void
  InitializeRegistrator();

private:
  std::string m_TheMovingImageFilename;
  std::string m_TheFixedImageFilename;
  std::string m_InitialDisplacementFieldFilename;
  std::string m_InitialCoefficientFilename;
  std::string m_InitialTransformFilename;
  std::string m_DisplacementBaseName;
  std::string m_WarpedImageName;

  // std::string m_ParameterFilename;
  bool m_ForceCoronalZeroOrigin;
  bool m_UseHistogramMatching;

  std::string m_OutNormalized;
  //    std::string m_OutDebug;
  std::string      m_OutputFilename;
  std::string      m_CheckerBoardFilename;
  std::string      m_DisplacementFieldOutputName;
  bool             m_AppendOutputFile;
  PatternArrayType m_CheckerBoardPattern;
  std::string      m_BOBFTargetMask;
  std::string      m_BOBFTemplateMask;
  IndexType        m_Seed;
  PixelType        m_Lower;
  PixelType        m_Upper;
  PixelType        m_DefaultPixelValue;
  SizeType         m_Radius; // for BOBF filter.
  SizeType         m_MedianFilterSize;
  std::string      m_FixedLandmarkFilename;
  std::string      m_MovingLandmarkFilename;
  // typename TDisplacementField::Pointer m_InitialDisplacementField;
  unsigned long       m_NumberOfHistogramLevels;
  unsigned long       m_NumberOfMatchPoints;
  unsigned short      m_NumberOfLevels;
  ShrinkFactorsType   m_TheMovingImageShrinkFactors;
  ShrinkFactorsType   m_TheFixedImageShrinkFactors;
  IterationsArrayType m_NumberOfIterations;
};
} // namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#  include "ThirionRegistration.hxx"
#endif
#endif
