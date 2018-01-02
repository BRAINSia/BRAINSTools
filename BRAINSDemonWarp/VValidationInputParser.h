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
#ifndef _VValidationInputParser_h
#define _VValidationInputParser_h

#include <string>
#include "itkObjectFactory.h"
#include "itkObject.h"
#include "itkFixedArray.h"
#include "itkArray.h"

namespace itk
{
/** \class VValidationInputParser
  *
  * This component parse an input parameter file for a simple
  * atlas based segmentation application.
  *
  * This class is activated by method Execute().
  *
  * Inputs:
  *  - altas image name
  *  - subject image name
  *  - the  parameter filename
  *
  *
  * Outputs:
  *  - pointer to the subject (fixed) image
  *  - pointer to the atlas (moving) image
  *  - the number of histogram levels to use
  *  - the number of histogram match points to use
  *  - the number of levels
  *  - the number of iterations at each level
  *  - the fixed image starting shrink factors
  *  - the moving image starting shrink factors
  *
  */

template <typename TImage>
class VValidationInputParser : public Object
{
public:

  /** Standard class typedefs. */
  typedef VValidationInputParser   Self;
  typedef Object                   Superclass;
  typedef SmartPointer<Self>       Pointer;
  typedef SmartPointer<const Self> ConstPointer;

  /** Run-time type information (and related methods). */
  itkTypeMacro(VValidationInputParser, Object);

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Image Type. */
  typedef TImage                      ImageType;
  typedef typename ImageType::Pointer ImagePointer;

  /** Image dimension enumeration. */
  itkStaticConstMacro(ImageDimension, unsigned int, TImage::ImageDimension);

  /** Displacement field value type. */
  typedef float FieldValueType;

  /** Displacement field pixel type. */
  typedef Vector<FieldValueType,
                 itkGetStaticConstMacro(ImageDimension)> FieldPixelType;

  /** Displacement field type. */
  typedef Image<FieldPixelType,
                itkGetStaticConstMacro(ImageDimension)> TDisplacementField;

  /** ShrinkFactors type. */
  typedef FixedArray<unsigned int,
                     itkGetStaticConstMacro(ImageDimension)> ShrinkFactorsType;

  /** IterationArray type. */
  typedef Array<unsigned int> IterationsArrayType;

  /** Set the atlas patient. */
  void SetTheMovingImageFilename(std::vector<std::string> & names)
  {
    m_TheMovingImageFilename = names;
  }

  /** Set the subject patient. */
  void SetTheFixedImageFilename(std::vector<std::string> & names)
  {
    m_TheFixedImageFilename = names;
  }

  /** Set the initial Displacement Field one of 3 ways. */
  itkSetStringMacro(InitialDisplacementFieldFilename);
  itkSetStringMacro(InitialCoefficientFilename);
  itkSetStringMacro(InitialTransformFilename);

  /** Set input parameter file name. */
  //            itkSetStringMacro( ParameterFilename );

  /** Parse the input file. */
  void Execute();

  /** Get pointer to the atlas image. */
  std::vector<ImagePointer> & GetTheMovingImages()
  {
    return m_TheMovingImages;
  }

  /** Get pointer to the subject image. */
  std::vector<ImagePointer> & GetTheFixedImages()
  {
    return m_TheFixedImages;
  }

  /**force Centered Image.*/
  itkSetMacro(ForceCoronalZeroOrigin, bool);
  itkGetConstMacro(ForceCoronalZeroOrigin, bool);

  /** Get pointer to the subject image. */
  itkGetModifiableObjectMacro(InitialDisplacementField, TDisplacementField);

  /** Get the number of histogram bins. */
  itkGetConstMacro(NumberOfHistogramLevels, unsigned long);
  itkSetMacro(NumberOfHistogramLevels, unsigned long);

  /** Get the number of match points. */
  itkGetConstMacro(NumberOfMatchPoints, unsigned long);
  itkSetMacro(NumberOfMatchPoints, unsigned long);

  /** Get the number of levels. */
  itkGetMacro(NumberOfLevels, unsigned short);
  itkSetMacro(NumberOfLevels, unsigned short);

  /**Set Debug mode*/
  itkSetMacro(OutDebug, bool);
  itkGetConstMacro(OutDebug, bool);

  /** Get the atlas image starting shrink factors. */
  itkGetConstReferenceMacro(TheMovingImageShrinkFactors, ShrinkFactorsType);
  void SetTheMovingImageShrinkFactors(const ShrinkFactorsType & shrinkfactors)
  {
    this->m_TheMovingImageShrinkFactors = shrinkfactors;
  }

  /** Get the subject image starting shrink factors. */
  itkGetConstReferenceMacro(TheFixedImageShrinkFactors, ShrinkFactorsType);
  void SetTheFixedImageShrinkFactors(const ShrinkFactorsType & shrinkfactors)
  {
    this->m_TheFixedImageShrinkFactors = shrinkfactors;
  }

  /** Get the number of iterations at each level. */
  itkGetConstReferenceMacro(NumberOfIterations, IterationsArrayType);
  void SetNumberOfIterations(const IterationsArrayType & iterations)
  {
    m_NumberOfIterations = iterations;
  }

protected:
  VValidationInputParser();
  ~VValidationInputParser() override
  {
  }

private:
  ITK_DISALLOW_COPY_AND_ASSIGN(VValidationInputParser);

  std::vector<std::string> m_TheMovingImageFilename;
  std::vector<std::string> m_TheFixedImageFilename;
  std::string              m_InitialDisplacementFieldFilename;
  std::string              m_InitialCoefficientFilename;
  std::string              m_InitialTransformFilename;
  std::string              m_ParameterFilename;

  std::vector<ImagePointer> m_TheMovingImages;
  std::vector<ImagePointer> m_TheFixedImages;

  bool m_ForceCoronalZeroOrigin;
  //  bool                          m_HistogramMatching;
  bool m_OutDebug;

  typename TDisplacementField::Pointer m_InitialDisplacementField;

  unsigned long       m_NumberOfHistogramLevels;
  unsigned long       m_NumberOfMatchPoints;
  unsigned short      m_NumberOfLevels;
  ShrinkFactorsType   m_TheMovingImageShrinkFactors;
  ShrinkFactorsType   m_TheFixedImageShrinkFactors;
  IterationsArrayType m_NumberOfIterations;
};
}   // namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "VValidationInputParser.hxx"
#endif

#endif
