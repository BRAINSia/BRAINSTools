/*=========================================================================

Program:   Insight Segmentation & Registration Toolkit
Module:    $RCSfile: itkHammerTissueAttributeVectorFromPartialVolumeImageFilter.h,v $
Language:  C++
Date:      $Date: 2009/01/13 20:19:20 $
Version:   $Revision: 1.4 $

Copyright (c) Insight Software Consortium. All rights reserved.
See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for
details.

This program is developed under NIH NCBC collaboration grant
R01 EB006733, "Development and Dissemination of Robust Brain MRI
Measurement Tools".

This software is distributed WITHOUT ANY WARRANTY; without even
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkHammerTissueAttributeVectorFromPartialVolumeImageFilter_h
#define __itkHammerTissueAttributeVectorFromPartialVolumeImageFilter_h

#include "itkImageToImageFilter.h"
#include "itkConstNeighborhoodIterator.h"
#include "itkHammerTissueAttributeVector.h"

namespace itk
{
/** \classHammerTissueAttributeVectorFromPartialVolumeImageFilter
  * \brief Computes the gradient of an image using directional derivatives.
  *
  * Computes the gradient of an image using directional derivatives.
  * The directional derivative at each pixel location is computed by
  * convolution with a first-order derivative operator.
  *
  * The second template parameter defines the value type used in the
  * derivative operator (defaults to float).  The third template
  * parameter defines the value type used for output image (defaults to
  * float).  The output image is defined as a covariant vector image
  * whose value type is specified as this third template parameter.
  *
  *
  * \sa Image
  * \sa Neighborhood
  * \sa NeighborhoodOperator
  * \sa NeighborhoodIterator
  *
  * \ingroup GradientFilters
  */
template <class TInputImage, class TOutputImage>
class HammerTissueAttributeVectorFromPartialVolumeImageFilter :
  public         ImageToImageFilter<TInputImage, TOutputImage>
{
public:
  /** Extract dimension from input image. */
  itkStaticConstMacro(InputImageDimension, unsigned int,
                      TInputImage::ImageDimension);
  itkStaticConstMacro(OutputImageDimension, unsigned int,
                      TInputImage::ImageDimension);

  /** Standard class typedefs. */
  typedef HammerTissueAttributeVectorFromPartialVolumeImageFilter Self;

  /** Convenient typedefs for simplifying declarations. */
  typedef TInputImage                          InputImageType;
  typedef typename InputImageType::Pointer     InputImagePointer;
  typedef TOutputImage                         OutputImageType;
  typedef typename OutputImageType::Pointer    OutputImagePointer;
  typedef typename InputImageType::IndexType   InputIndexType;
  typedef typename InputImageType::SpacingType InputSpacingType;
  typedef typename InputImageType::RegionType  InputRegionType;

  /** types for neighborhood iterator */
  typedef ConstNeighborhoodIterator<InputImageType>     NeighborhoodIteratorType;
  typedef typename NeighborhoodIteratorType::OffsetType NeighborOffsetType;

  /** Standard class typedefs. */
  typedef ImageToImageFilter<InputImageType, OutputImageType> Superclass;
  typedef SmartPointer<Self>                                  Pointer;
  typedef SmartPointer<const Self>                            ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(HammerTissueAttributeVectorFromPartialVolumeImageFilter, ImageToImageFilter);

  /** Image typedef support. */
  typedef typename InputImageType::PixelType   InputPixelType;
  typedef typename OutputImageType::PixelType  OutputValueType;
  typedef typename OutputImageType::RegionType OutputImageRegionType;

  /** Use the image spacing information in calculations. Use this option if you
    *  want derivatives in physical space. Default is UseImageSpacingOn. */
  void SetUseImageSpacingOn()
  {
    this->SetUseImageSpacing(true);
  }

  /** Ignore the image spacing. Use this option if you want derivatives in
    isotropic pixel space.  Default is UseImageSpacingOn. */
  void SetUseImageSpacingOff()
  {
    this->SetUseImageSpacing(false);
  }

  // set input prior images
  //
  // 0. inputGMVolume  : Grey Matter Partial volume/Posterior image
  // 1. inputWMVolume  : White Matter Partial volume/Posteior image
  // 3. inputCSFVolume : CSF Partial volume/Posterior image

  void SetGMVolume( const TInputImage * inputGMVolume)
  {
    this->SetNthInput(0, const_cast<TInputImage *>( inputGMVolume ) );
  }

  void SetWMVolume( const TInputImage * inputWMVolume)
  {
    this->SetNthInput(1, const_cast<TInputImage *>( inputWMVolume ) );
  }

  void SetCSFVolume( const TInputImage * inputCSFVolume)
  {
    this->SetNthInput(2, const_cast<TInputImage *>( inputCSFVolume ) );
  }

  /** HammerTissueAttributeVectorFromPartialVolumeImageFilter needs a larger input requested region than
    * the output requested region.  As such, GradientImageFilter needs
    * to provide an implementation for GenerateInputRequestedRegion()
    * in order to inform the pipeline execution model.
    *
    * \sa ImageToImageFilter::GenerateInputRequestedRegion() */
  virtual void GenerateInputRequestedRegion()
  throw ( InvalidRequestedRegionError );

  /** Set/Get whether or not the filter will use the spacing of the input
    image in its calculations */
  itkSetMacro(UseImageSpacing, bool);
  itkGetMacro(UseImageSpacing, bool);

  /** Set/Get macroes for class member varibles */
  itkSetMacro(Scale, float);
  itkGetMacro(Scale, float);

  itkSetMacro(Strength, unsigned char);
  itkGetMacro(Strength, unsigned char);

  itkSetMacro(GMValue, unsigned char);
  itkGetMacro(GMValue, unsigned char);

  itkSetMacro(WMValue, unsigned char);
  itkGetMacro(WMValue, unsigned char);

  itkSetMacro(CSFValue, unsigned char);
  itkGetMacro(CSFValue, unsigned char);

  itkSetMacro(VNValue, unsigned char);
  itkGetMacro(VNValue, unsigned char);

  itkSetMacro(BGValue, unsigned char);
  itkGetMacro(BGValue, unsigned char);

  /** The UseImageDirection flag determines whether image derivatives are
    * computed with respect to the image grid or with respect to the physical
    * space. When this flag is ON the derivatives are computed with respect to
    * the coodinate system of physical space. The difference is whether we take
    * into account the image Direction or not. The flag ON will take into
    * account the image direction and will result in an extra matrix
    * multiplication compared to the amount of computation performed when the
    * flag is OFF.
    * The default value of this flag is the same as the CMAKE option
    * ITK_IMAGE_BEHAVES_AS_ORIENTED_IMAGE (i.e ON by default when ITK_IMAGE_BEHAVES_AS_ORIENTED_IMAGE is ON,
    * and  OFF by default when ITK_IMAGE_BEHAVES_AS_ORIENTED_IMAGE is OFF).*/
  itkSetMacro(UseImageDirection, bool);
  itkGetMacro(UseImageDirection, bool);
  itkBooleanMacro(UseImageDirection);
protected:
  HammerTissueAttributeVectorFromPartialVolumeImageFilter();
  virtual ~HammerTissueAttributeVectorFromPartialVolumeImageFilter();
  void PrintSelf(std::ostream & os, Indent indent) const;

  /** HammerTissueAttributeVectorFromPartialVolumeImageFilter can be implemented as a multithreaded filter.
  * Therefore, this implementation provides a ThreadedGenerateData()
  * routine which is called for each processing thread. The output
  * image data is allocated automatically by the superclass prior to
  * calling ThreadedGenerateData().  ThreadedGenerateData can only
  * write to the portion of the output image specified by the
    * parameter "outputRegionForThread"
    *
    * \sa ImageToImageFilter::ThreadedGenerateData(),
    *     ImageToImageFilter::GenerateData() */
  virtual void GenerateData();

private:
  HammerTissueAttributeVectorFromPartialVolumeImageFilter(const Self &);   // purposely not
  // implemented
  void operator=(const Self &);                           // purposely not

  // implemented

  void CreateN1Neighbor();

  void CreateFeatureNeighbor(int Radius);

  bool m_UseImageSpacing;

  // flag to take or not the image direction into account
  // when computing the derivatives.
  bool m_UseImageDirection;

  // scale for computing the geometric features
  float m_Scale;

  // edge "stength" for computing edges
  unsigned char m_Strength;

  // labels for Gray Matter, White Matter, Cortical CSF, Ventricular
  // CSF, and Background
  unsigned char m_GMValue;
  unsigned char m_WMValue;
  unsigned char m_CSFValue;
  unsigned char m_VNValue;
  unsigned char m_BGValue;

  // const values for edge type
  static const unsigned char m_WMGMEDGE = 255;
  static const unsigned char m_WMCSFEDGE = 180;
  static const unsigned char m_GMCSFEDGE = 100;

  std::vector<NeighborOffsetType> m_FeatureNeighborhood;
  std::vector<InputIndexType>     m_N1Neighborhood;
  // indices in the spherical neighborhood
  std::vector<NeighborOffsetType> m_OffsetInSphericalNeighborhood;
};
}   // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkHammerTissueAttributeVectorFromPartialVolumeImageFilter.hxx"
#endif

#endif
