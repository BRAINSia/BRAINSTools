/*
 *  itkResampleInPlaceImageFilter.h
 *
 *
 *  Created by Wei Lu on 10/14/10.
 *
 */

#ifndef __itkResampleInPlaceImageFilter_h
#define __itkResampleInPlaceImageFilter_h

#include "itkImageToImageFilter.h"
#include "itkVersorRigid3DTransform.h"

namespace itk
{
/** \class ResampleInPlaceImageFilter
 * \brief Resample an image in place.
 *
 * The current ITK resample image filter will generate a physical memory-
 * modified version of the input image if the input transform is not identical. The
 * abuse use of the filter can be cumbersome in a situation that the image is very
 * large, and there are lots of transform to be superimposed for the input image, and
 * we even don't care about those intermediate transformed images.
 *
 * If all the transforms are rigid, a far superior way to achieve a similar result
 * while minimizing the accumulative resampling errors as well as eliminating the expense
 * on accessing the physical memory of the image is to compose all the
 * transforms before hand if it is possible, and then we only need to resample the
 * input image with the final composed transform once.
 *
 * Here we present a more compact alternative – all information is stored in the header
 * of the image and there is no need to maintain the final transform any longer. ITK
 * image class has innate support for doing this.
 *
 * \param RigidTransform -- Currently must be a VersorRigid3D
 * \param InputImage -- The image to be duplicated and modified to incorporate the
 * rigid transform
 * \return -- An image with the same voxels values as the input, but with differnt
 * physical space representation affected by the rigid transform.
 *
 * \ingroup GeometricTransforms
 */
template <class TInputImage, class TOutputImage>
class ITK_EXPORT ResampleInPlaceImageFilter :
  public         ImageToImageFilter<TInputImage, TOutputImage>
{
public:
  /** Standard class typedefs */
  typedef ResampleInPlaceImageFilter                    Self;
  typedef ImageToImageFilter<TInputImage, TOutputImage> Superclass;
  typedef SmartPointer<Self>                            Pointer;
  typedef SmartPointer<const Self>                      ConstPointer;

  /** Method for creation through the object factory */
  itkNewMacro( Self );

  /** Run-time type information (and related methods) */
  itkTypeMacro( ResampleInPlaceImageFilter, ImageToImageFilter );

  /** input/output image typedefs */
  typedef TInputImage                         InputImageType;
  typedef typename InputImageType::Pointer    InputImagePointer;
  typedef typename InputImageType::RegionType InputImageRegionType;
  typedef typename InputImageType::PixelType  InputImagePixelType;
  typedef typename InputImageType::PointType  InputImagePointType;

  typedef TOutputImage                         OutputImageType;
  typedef typename OutputImageType::Pointer    OutputImagePointer;
  typedef typename OutputImageType::RegionType OutputImageRegionType;
  typedef typename OutputImageType::PixelType  OutputImagePixelType;

  /** ImageDimension constants */
  itkStaticConstMacro(InputImageDimension, unsigned int,
                      TInputImage::ImageDimension);
  itkStaticConstMacro(OutputImageDimension, unsigned int,
                      TOutputImage::ImageDimension);

#ifdef ITK_USE_CONCEPT_CHECKING
  itkConceptMacro(SameDimensionCheck,
                  (Concept::SameDimension<InputImageDimension, OutputImageDimension> ) );
  itkConceptMacro(InputConvertibleToOutputCheck,
                  (Concept::Convertible<InputImagePixelType, OutputImagePixelType> ) );
#endif

  /** Transform typedef */
  typedef VersorRigid3DTransform<double>            RigidTransformType;
  typedef typename RigidTransformType::ConstPointer RigidTransformConstPointer;

  /** Set/Get rigid transform. The default is an identity transform */
  itkSetConstObjectMacro( RigidTransform, RigidTransformType );
  itkGetConstObjectMacro( RigidTransform, RigidTransformType );

  /** Set/Get required input image. (A wrapper to this->Set/GetInput()) */
  void SetInputImage( const InputImageType * image );

  const InputImageType * GetInputImage() const;

protected:
  ResampleInPlaceImageFilter();
  ~ResampleInPlaceImageFilter()
  {
  };

  void GenerateData();

  void PrintSelf( std::ostream& os, Indent indent ) const;

private:
  ResampleInPlaceImageFilter( const Self & ); // purposely not implemented
  void operator=( const Self & );             // purposely not implemented

  OutputImagePointer         m_OutputImage;
  RigidTransformConstPointer m_RigidTransform;
};
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkResampleInPlaceImageFilter.hxx"
#endif

#endif
