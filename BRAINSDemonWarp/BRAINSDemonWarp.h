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
#ifndef _BRAINSDemonWarp_h
#define _BRAINSDemonWarp_h

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

template <typename TImage,
          typename TRealImage, typename TOutputImage
          >
class BRAINSDemonWarp : public ApplicationBase<
    ValidationInputParser<TImage>,
    DemonsPreprocessor<TImage, TRealImage>,
    DemonsRegistrator<TRealImage, TOutputImage,
                      typename TRealImage::PixelType>
    >
{
public:

  /** Standard class typedefs. */
  typedef BRAINSDemonWarp Self;
  typedef ApplicationBase<ValidationInputParser<TImage>,
                          DemonsPreprocessor<TImage, TRealImage>,
                          DemonsRegistrator<TRealImage, TRealImage,
                                            typename TRealImage::PixelType>
                          > Superclass;
  typedef SmartPointer<Self>       Pointer;
  typedef SmartPointer<const Self> ConstPointer;

  /** Displacement field pixel type. */
  typedef float                     FieldValueType;
  typedef Vector<FieldValueType, 3> FieldPixelType;
  typedef Image<FieldPixelType, 3>  TDisplacementField;

  /** Run-time type information (and related methods). */
  itkTypeMacro(BRAINSDemonWarp, ApplicationBase);

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Image types. */
  typedef TImage     ImageType;
  typedef TRealImage RealImageType;

  /** Image dimension. */
  itkStaticConstMacro(ImageDimension, unsigned int,
                      TImage::ImageDimension);

  /** Type to hold the number of checker boxes per dimension */
  typedef FixedArray<unsigned int, TImage::ImageDimension> PatternArrayType;

  typedef typename ImageType::PixelType PixelType;
  typedef typename ImageType::IndexType IndexType;
  typedef typename ImageType::SizeType  SizeType;

  /** ShrinkFactors type. */
  typedef FixedArray<unsigned int, TImage::ImageDimension>
    ShrinkFactorsType;

  /** IterationArray type. */
  typedef Array<unsigned int> IterationsArrayType;

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

  /**Set checker board Image filename*/
  itkSetStringMacro(CheckerBoardFilename);
  itkGetStringMacro(CheckerBoardFilename);

  /**Set Displacement field output filename*/
  itkSetStringMacro(DisplacementFieldOutputName);
  itkGetStringMacro(DisplacementFieldOutputName);

  /** Set Checker pattern */
  itkSetMacro(CheckerBoardPattern, PatternArrayType);
  itkGetConstReferenceMacro(CheckerBoardPattern, PatternArrayType);

  /** Set append output file boolean. */
  itkSetMacro(AppendOutputFile, bool);
  itkGetMacro(AppendOutputFile, bool);
  itkBooleanMacro(AppendOutputFile);

  /*BOBF macros
   *Set Target Mask filename*/
  itkSetStringMacro(FixedBinaryVolume);
  itkGetStringMacro(FixedBinaryVolume);

  /**Set Template Mask filename*/
  itkSetStringMacro(MovingBinaryVolume);
  itkGetStringMacro(MovingBinaryVolume);

  /**Force Centered Image.*/
  itkSetMacro(ForceCoronalZeroOrigin, bool);
  itkGetConstMacro(ForceCoronalZeroOrigin, bool);

  /**Output Normalized Image.*/
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

  itkSetMacro(MedianFilterSize,  SizeType);
  itkGetMacro(MedianFilterSize,  SizeType);

  /** Set the initial deformation field to prime registration */
  //    itkSetObjectMacro(InitialDisplacementField,TDisplacementField);
  /** Set the Input Landmark Filename*/
  itkSetStringMacro(FixedLandmarkFilename);
  itkGetStringMacro(FixedLandmarkFilename);
  itkSetStringMacro(MovingLandmarkFilename);
  itkGetStringMacro(MovingLandmarkFilename);
  /*
    * itkSetStringMacro (FixedMaskFilename);
    * itkGetStringMacro (FixedMaskFilename);
    * itkSetStringMacro (MovingMaskFilename);
    * itkGetStringMacro (MovingMaskFilename);
    */
  /**Set histogram matching*/
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

  /** Get the interpolation Mode. */
  itkGetMacro(InterpolationMode, std::string);
  itkSetMacro(InterpolationMode, std::string);

  /** Get the atlas image starting shrink factors. */
  itkGetConstReferenceMacro(TheMovingImageShrinkFactors, ShrinkFactorsType);
  void SetTheMovingImageShrinkFactors(const ShrinkFactorsType & shrinkfactors)
  {
    this->m_TheMovingImageShrinkFactors = shrinkfactors;
    this->Modified();
  }

  /** Get the subject image starting shrink factors. */
  itkGetConstReferenceMacro(TheFixedImageShrinkFactors, ShrinkFactorsType);
  void SetTheFixedImageShrinkFactors(const ShrinkFactorsType & shrinkfactors)
  {
    this->m_TheFixedImageShrinkFactors = shrinkfactors;
    this->Modified();
  }

  /** Get the number of iterations at each level. */
  itkGetConstReferenceMacro(NumberOfIterations, IterationsArrayType);
  void SetNumberOfIterations(const IterationsArrayType & iterations)
  {
    m_NumberOfIterations = iterations;
    this->Modified();
  }

  typedef itk::PDEDeformableRegistrationFilter<RealImageType, RealImageType,
                                               TDisplacementField>
    BaseRegistrationFilterType;
  void SetRegistrationFilter(
    BaseRegistrationFilterType *filter)
  {
    this->m_Registrator->SetRegistrationFilter(filter);
    this->Modified();
  }

protected:

  BRAINSDemonWarp();
  ~BRAINSDemonWarp() override
  {
  }

  /** Initialize the input parser. */
  void InitializeParser() override;

  /*** Initialize the preprocessor */
  void InitializePreprocessor() override;

  /*** Initialize the registrator  */
  void InitializeRegistrator() override;

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
  std::string      m_FixedBinaryVolume;
  std::string      m_MovingBinaryVolume;
  IndexType        m_Seed;
  PixelType        m_Lower;
  PixelType        m_Upper;
  PixelType        m_DefaultPixelValue;
  SizeType         m_Radius;   // for BOBF filter.
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
  std::string         m_InterpolationMode;
};
}          // namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "BRAINSDemonWarp.hxx"
#endif
#endif
