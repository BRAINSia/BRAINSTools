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
#ifndef __DemonsPreprocessor_h
#define __DemonsPreprocessor_h

#include "itkObject.h"

namespace itk
{
/** \class DemonsPreprocessor
  *
  * This component pre-processes the moving and fixed image before
  * registration.
  * If the fixed image dimensions are different from the moving image it will
  *    * resample the moving image to match the fixed image dimensions.
  * Histogram matching is done to solve the intensity mismatch problem.
  *
  * The preprocessor also called the skull-stripping filter itkBOBF
  * if an atlas and subject whole brain masks are specified.
  *
  * The preprocessing is activatived by method Execute().
  *
  * Inputs:
  *    - pointer to original fixed image
  *    - pointer original moving image
  *    - number of histogram levels
  *    - number of match points
  *
  * Outputs:
  *    - pointer to processed fixed image
  *    - pointer to processed moving image
  *    - the minimum value of original fixed image
  *    - the minimum value of original moving image
  *
  */
template <typename TInputImage, typename TOutputImage>
class DemonsPreprocessor : public Object
{
public:

  /** Standard class typedefs. */
  typedef DemonsPreprocessor       Self;
  typedef Object                   Superclass;
  typedef SmartPointer<Self>       Pointer;
  typedef SmartPointer<const Self> ConstPointer;

  /** Run-time type information (and related methods). */
  itkTypeMacro(DemonsPreprocessor, Object);

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Input Image Type. */
  typedef TInputImage InputImageType;
  /** Output Image Type. */
  typedef TOutputImage OutputImageType;

  /** Input image pixel type. */
  typedef typename InputImageType::PixelType  InputPixelType;
  typedef typename OutputImageType::PixelType PixelType;
  typedef typename OutputImageType::IndexType IndexType;
  typedef typename OutputImageType::SizeType  SizeType;

  /** Image dimension enumeration. */
  itkStaticConstMacro(ImageDimension, unsigned int, TInputImage::ImageDimension);

  /** Set the input fixed image. */
  itkSetObjectMacro(InputFixedImage, InputImageType);

  /** Set the input moving image. */
  itkSetObjectMacro(InputMovingImage, InputImageType);

  /** Displacement field value type. */
  typedef float FieldValueType;

  /** Displacement field pixel type. */
  typedef Vector<FieldValueType,
                 itkGetStaticConstMacro(ImageDimension)> FieldPixelType;

  /** Displacement field type. */
  typedef Image<FieldPixelType,
                itkGetStaticConstMacro(ImageDimension)> TDisplacementField;

  /** Set the initial Displacement Field. */
  itkSetObjectMacro(InitialDisplacementField, TDisplacementField);
  itkGetModifiableObjectMacro(InitialDisplacementField, TDisplacementField);

  /** Set the number of histogram levels to use. */
  itkSetMacro(NumberOfHistogramLevels, unsigned long);

  /** Set the number of match points to use. */
  itkSetMacro(NumberOfMatchPoints, unsigned long);

  /** Method to execute the preprocessing. */
  virtual void Execute();

  /** Get the output fixed/moving image. */
  itkGetModifiableObjectMacro(OutputFixedImage, OutputImageType);
  itkGetModifiableObjectMacro(OutputMovingImage, OutputImageType);

  /** Get the output fixed/moving image. */
  itkGetModifiableObjectMacro(UnNormalizedMovingImage, OutputImageType);
  itkGetModifiableObjectMacro(UnNormalizedFixedImage, OutputImageType);

  /** Get minimum value of original fixed image. */
  itkGetMacro(FixedImageMinimum, InputPixelType);

  /** Get minimum value of original moving image. */
  itkGetMacro(MovingImageMinimum, InputPixelType);

  /* BOBF macros
    * Set Target Mask filename */
  itkSetStringMacro(FixedBinaryVolume);
  itkGetStringMacro(FixedBinaryVolume);

  /** Set Template Mask filename */
  itkSetStringMacro(MovingBinaryVolume);
  itkGetStringMacro(MovingBinaryVolume);

  /** Set/Get the lower threshold. The default is 0. */
  itkSetMacro(Lower, PixelType);
  itkGetMacro(Lower, PixelType);

  /** Set/Get the upper threshold. The default is 70 */
  itkSetMacro(Upper, PixelType);
  itkGetMacro(Upper, PixelType);

  itkSetMacro(DefaultPixelValue,  PixelType);
  itkGetMacro(DefaultPixelValue,  PixelType);

  itkSetMacro(MedianFilterSize,  SizeType);
  itkGetMacro(MedianFilterSize,  SizeType);

  /** Set the radius of the neighborhood used for a mask. */
  itkSetMacro(Radius, SizeType);
  /** Get the radius of the neighborhood used to compute the median */
  itkGetConstReferenceMacro(Radius, SizeType);

  /** Set the Seed of the neighborhood used for a mask. */
  itkSetMacro(Seed, IndexType);
  /** Get the radius of the neighborhood used to compute the median */
  itkGetConstReferenceMacro(Seed, IndexType);

  /**Set Debug mode*/
  itkSetMacro(OutDebug, bool);
  itkGetConstMacro(OutDebug, bool);

  /**Set histogram matching*/
  itkSetMacro(UseHistogramMatching, bool);
  itkGetConstMacro(UseHistogramMatching, bool);
protected:
  DemonsPreprocessor();
  ~DemonsPreprocessor()
  {
  }

private:
  ITK_DISALLOW_COPY_AND_ASSIGN(DemonsPreprocessor);

  typename InputImageType::Pointer m_InputFixedImage;
  typename InputImageType::Pointer m_InputMovingImage;
  typename OutputImageType::Pointer m_OutputFixedImage;
  typename OutputImageType::Pointer m_OutputMovingImage;
  typename OutputImageType::Pointer m_UnNormalizedMovingImage;
  typename OutputImageType::Pointer m_UnNormalizedFixedImage;
  typename TDisplacementField::Pointer m_InitialDisplacementField;

  unsigned long m_NumberOfHistogramLevels;
  unsigned long m_NumberOfMatchPoints;

  InputPixelType m_FixedImageMinimum;
  InputPixelType m_MovingImageMinimum;

  std::string m_FixedBinaryVolume;
  std::string m_MovingBinaryVolume;
  IndexType   m_Seed;
  PixelType   m_Lower;
  PixelType   m_Upper;
  PixelType   m_DefaultPixelValue;
  SizeType    m_Radius;
  bool        m_OutDebug;
  SizeType    m_MedianFilterSize;
  typedef typename OutputImageType::Pointer OutputImagePointer;
  typedef typename InputImageType::Pointer  InputImagePointer;

  bool m_UseHistogramMatching;

  /*MakeBOBF function takes in a brain image and a whole brain mask and strips
    * the skull of the image.*/
  OutputImagePointer MakeBOBFImage(OutputImagePointer input, std::string MaskName);
};
}   // namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "DemonsPreprocessor.hxx"
#endif

#endif // _DemonsPreprocessor_h
