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
#ifndef _IccdefPreprocessor_h
#define _IccdefPreprocessor_h

#include "itkObject.h"
#include "itkAffineTransform.h"
#include "itkScaleVersor3DTransform.h"

#include "itkVersorRigid3DTransform.h"
#include "itkScaleSkewVersor3DTransform.h"
#include "itkBSplineTransform.h"

namespace itk
{
/** \class IccdefPreprocessor
 *
 * This component pre-processes the moving and fixed image before
 * registration.
 * If the fixed image dimensions are different from the moving image it will     * resample the moving image to match the fixed image dimensions.
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
class IccdefPreprocessor : public Object
{
public:

  /** Standard class typedefs. */
  typedef IccdefPreprocessor       Self;
  typedef Object                   Superclass;
  typedef SmartPointer<Self>       Pointer;
  typedef SmartPointer<const Self> ConstPointer;

  /** Run-time type information (and related methods). */
  itkTypeMacro(IccdefPreprocessor, Object);

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

  /** Set the atlas patient ID. */
  itkSetStringMacro(TheMovingImageFilename);
  itkGetStringMacro(TheMovingImageFilename);

  /** Set the subject patient ID. */
  itkSetStringMacro(TheFixedImageFilename);
  itkGetStringMacro(TheFixedImageFilename);

  /** Deformation field value type. */
  typedef float FieldValueType;

  typedef itk::Image<unsigned char, itkGetStaticConstMacro(ImageDimension)> TImageType;

  /** Deformation field pixel type. */
  typedef Vector<FieldValueType,
                 itkGetStaticConstMacro(ImageDimension)> FieldPixelType;

  /** Deformation field type. */
  typedef Image<FieldPixelType,
                itkGetStaticConstMacro(ImageDimension)> TDisplacementField;

  /** Set the initial Deformation Field. */
  itkSetObjectMacro( InitialDisplacementField, TDisplacementField );
  itkGetModifiableObjectMacro( InitialDisplacementField, TDisplacementField );

  /** Set the number of histogram levels to use. */
  itkSetMacro( NumberOfHistogramLevels, unsigned long );

  /** Set the number of match points to use. */
  itkSetMacro( NumberOfMatchPoints, unsigned long );

  /** Method to execute the preprocessing. */
  virtual void Execute();

  /** Get the output fixed image. */
  itkGetModifiableObjectMacro( OutputFixedImage, OutputImageType );
  itkGetModifiableObjectMacro( OutputMovingImage, OutputImageType );

  /** Get the output moving image. */
  itkGetModifiableObjectMacro( UnNormalizedMovingImage, OutputImageType );
  itkGetModifiableObjectMacro( UnNormalizedFixedImage, OutputImageType );

  /** Get minimum value of original fixed image. */
  itkGetMacro( FixedImageMinimum, InputPixelType );

  /** Get minimum value of original moving image. */
  itkGetMacro( MovingImageMinimum, InputPixelType );

  /** Transform Types. */
  typedef VersorRigid3DTransform<double>                                  VersorRigid3DTransformType;
  typedef ScaleSkewVersor3DTransform<double>                              ScaleSkewVersor3DTransformType;
  typedef AffineTransform<double, itkGetStaticConstMacro(ImageDimension)> AffineTransformType;

  /** Image dimension enumeration. */
  itkStaticConstMacro(SplineOrder, unsigned int, 3);

  typedef double CoordinateRepType;
  typedef typename itk::BSplineTransform<
      CoordinateRepType,
      itkGetStaticConstMacro(ImageDimension),
      itkGetStaticConstMacro(SplineOrder)> BSplineTransformType;

  /*BOBF macros*/
  /**Set Target Mask filename*/
  itkSetStringMacro(FixedBinaryVolume );
  itkGetStringMacro( FixedBinaryVolume );

  /**Set Template Mask filename*/
  itkSetStringMacro(MovingBinaryVolume );
  itkGetStringMacro( MovingBinaryVolume );


  itkSetMacro(DefaultPixelValue,  PixelType);
  itkGetMacro(DefaultPixelValue,  PixelType);

  itkSetMacro(MedianFilterSize,  SizeType);
  itkGetMacro(MedianFilterSize,  SizeType);

  itkSetStringMacro( InitialTransformFilename );

  /**Set Debug mode*/
  itkSetMacro(OutDebug, bool);
  itkGetConstMacro(OutDebug, bool);

  /**Set histogram matching*/
  itkSetMacro(UseHistogramMatching, bool);
  itkGetConstMacro(UseHistogramMatching, bool);
protected:
  IccdefPreprocessor();
  ~IccdefPreprocessor()
  {
  };
private:
  ITK_DISALLOW_COPY_AND_ASSIGN(IccdefPreprocessor);

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
//  IndexType   m_Seed;
//  PixelType   m_Lower;
//  PixelType   m_Upper;
//  SizeType    m_Radius;
  PixelType m_DefaultPixelValue;
  bool      m_OutDebug;
  SizeType  m_MedianFilterSize;
  typedef typename OutputImageType::Pointer OutputImagePointer;
  typedef typename InputImageType::Pointer  InputImagePointer;

  bool        m_UseHistogramMatching;
  std::string m_TheMovingImageFilename;
  std::string m_TheFixedImageFilename;
  std::string m_InitialTransformFilename;

  /*MakeBOBF function takes in a brain image and a whole brain mask and strips
    the skull of the image.*/
  InputImagePointer MakeBOBFImage( InputImagePointer input, std::string MaskName );
};
}   // namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "IccdefPreprocessor.hxx"
#endif

#endif // _IccdefPreprocessor_h
