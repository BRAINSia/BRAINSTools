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
/**
  * \defgroup CF Create Field
  * \ingroup Reg
  */
#ifndef __CreateField_h
#define __CreateField_h

#include "itkObjectFactory.h"
#include "itkObject.h"
#include "itkFixedArray.h"
#include "itkArray.h"
#include "itkVector.h"
#include "itkImage.h"
#include "itkMultiResolutionPDEDeformableRegistration.h"
#include "itkRecursiveMultiResolutionPyramidImageFilter.h"

namespace itk
{
template <typename TImage, typename T2Image>
class CreateField : public Object
{
public:
  typedef CreateField              Self;
  typedef Object                   Superclass;
  typedef SmartPointer<Self>       Pointer;
  typedef SmartPointer<const Self> ConstPointer;

  itkTypeMacro(MIMApplication, Object);

  itkNewMacro(Self);

  itkSetStringMacro(Image1Filename);
  itkGetStringMacro(Image1Filename);
  itkSetStringMacro(Image2Filename);
  itkGetStringMacro(Image2Filename);
  itkSetStringMacro(ParameterFilename);

  typedef TImage                      ImageType;
  typedef typename ImageType::Pointer ImagePointer;
  static constexpr unsigned int ImageDimension = TImage::ImageDimension;
  typedef Array<unsigned int> IterationsArrayType;
  itkGetConstObjectMacro(ImageOne, ImageType);
  itkGetConstObjectMacro(ImageTwo, ImageType);
  itkSetObjectMacro(ImageOne, ImageType);
  itkSetObjectMacro(ImageTwo, ImageType);
  itkSetMacro(NumberOfHistogramLevels, unsigned long);
  itkGetMacro(NumberOfHistogramLevels, unsigned long);

  itkGetMacro(NumberOfMatchPoints, unsigned long);
  itkSetMacro(NumberOfMatchPoints, unsigned long);

  itkGetMacro(NumberOfLevels, unsigned short);
  typedef FixedArray<unsigned int,
                     itkGetStaticConstMacro(ImageDimension)> ShrinkFactorsType;
  itkGetMacro(Image1ShrinkFactors, ShrinkFactorsType);
  itkSetMacro(Image1ShrinkFactors, ShrinkFactorsType);
  itkGetMacro(Image2ShrinkFactors, ShrinkFactorsType);
  itkSetMacro(Image2ShrinkFactors, ShrinkFactorsType);
  itkGetConstReferenceMacro(NumberOfIterations, IterationsArrayType);

  typedef TImage                             InputImageType;
  typedef typename InputImageType::PixelType InputPixelType;
  typedef T2Image                            OutputImageType;
  itkGetConstObjectMacro(FixedImage, OutputImageType);
  itkGetConstObjectMacro(MovingImage, OutputImageType);
  itkGetMacro(FixedImageMinimum, InputPixelType);
  itkGetMacro(MovingImageMinimum, InputPixelType);

  typedef TImage
    FixedImageType;
  typedef T2Image
    MovingImageType;
  typedef Vector<float,
                 itkGetStaticConstMacro(ImageDimension)> FieldPixelType;
  typedef Image<FieldPixelType,
                itkGetStaticConstMacro(ImageDimension)> TDisplacementField;
  typedef RecursiveMultiResolutionPyramidImageFilter<FixedImageType,
                                                     FixedImageType>                         FixedImagePyramidType;
  typedef RecursiveMultiResolutionPyramidImageFilter<MovingImageType,
                                                     MovingImageType>                        MovingImagePyramidType;
  typedef MultiResolutionPDEDeformableRegistration<FixedImageType,
                                                   MovingImageType, TDisplacementField>     RegistrationType;

  typedef Array<unsigned int>
    UnsignedIntArray;
  itkSetClampMacro( NumberOfLevels, unsigned short, 1,
                    NumericTraits<unsigned short>::max() );
  itkSetMacro(NumberOfIterations, UnsignedIntArray);
  itkGetConstObjectMacro(DisplacementField, TDisplacementField);
  void StartNewLevel();

  void Execute();

  void ReleaseDataFlagOn();

protected:
  CreateField();
  virtual ~CreateField();
private:
  typename ImageType::Pointer m_ImageOne;
  typename ImageType::Pointer m_ImageTwo;
  std::string m_Image1Filename;
  std::string m_Image2Filename;
  std::string m_ParameterFilename;

  unsigned long     m_NumberOfHistogramLevels;
  unsigned long     m_NumberOfMatchPoints;
  unsigned short    m_NumberOfLevels;
  ShrinkFactorsType m_Image1ShrinkFactors;
  ShrinkFactorsType m_Image2ShrinkFactors;
  UnsignedIntArray  m_NumberOfIterations;

  InputPixelType m_FixedImageMinimum;
  InputPixelType m_MovingImageMinimum;

  typename OutputImageType::Pointer m_FixedImage;
  typename OutputImageType::Pointer m_MovingImage;
  typename FixedImagePyramidType::Pointer m_FixedImagePyramid;
  typename MovingImagePyramidType::Pointer m_MovingImagePyramid;
  typename TDisplacementField::Pointer m_DisplacementField;
  unsigned long m_Tag;
  typename RegistrationType::Pointer m_Registration;

  typedef typename OutputImageType::Pointer OutputImagePointer;
  void NormalizeImage(InputImageType *input, OutputImagePointer & output, InputPixelType & min);
};
}

#ifndef ITK_MANUAL_INSTANTIATION
#include "CreateField.hxx"
#endif

#endif
