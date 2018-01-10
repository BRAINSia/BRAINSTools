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
/*=========================================================================
 *
 *  Program:   Insight Segmentation & Registration Toolkit
 *  Module:    $RCSfile: itkBRAINSROIAutoImageFilter.h,v $
 *  Language:  C++
 *  Date:      $Date: 2008-10-16 19:33:40 $
 *  Version:   $Revision: 1.7 $
 *
 *  Copyright (c) Insight Software Consortium. All rights reserved.
 *  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.
 *
 *  This software is distributed WITHOUT ANY WARRANTY; without even
 *  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 *  PURPOSE.  See the above copyright notices for more information.
 *
 *  =========================================================================*/
#ifndef __itkBRAINSROIAutoImageFilter_h
#define __itkBRAINSROIAutoImageFilter_h

// First make sure that the configuration is available.
// This line can be removed once the optimized versions
// gets integrated into the main directories.
#include "itkConfigure.h"
#include "itkImageToImageFilter.h"
#include "itkImage.h"

#include "BRAINSTypes.h"
#include "itkLargestForegroundFilledMaskImageFilter.h"
#include "itkCastImageFilter.h"

namespace itk
{
/** \class BRAINSROIAutoImageFilter
  * \brief This is a class to help with identifying common tissue
  * Regions in an image.
  *
  * \sa Image
  * \sa Neighborhood
  *
  * \ingroup IntensityImageFilters
  */
template <class TInputImage, class TOutputImage>
class BRAINSROIAutoImageFilter :
  public         ImageToImageFilter<TInputImage, TOutputImage>
{
public:
  /** Extract dimension from input and output image. */
  itkStaticConstMacro(InputImageDimension, unsigned int,
                      TInputImage::ImageDimension);
  itkStaticConstMacro(OutputImageDimension, unsigned int,
                      TOutputImage::ImageDimension);

  /** Convenient typedefs for simplifying declarations. */
  typedef TInputImage  InputImageType;
  typedef TOutputImage OutputImageType;

  /** Standard class typedefs. */
  typedef BRAINSROIAutoImageFilter                            Self;
  typedef ImageToImageFilter<InputImageType, OutputImageType> Superclass;
  typedef SmartPointer<Self>                                  Pointer;
  typedef SmartPointer<const Self>                            ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(BRAINSROIAutoImageFilter, ImageToImageFilter);

  /** Image typedef support. */
  typedef typename InputImageType::PixelType  InputPixelType;
  typedef typename OutputImageType::PixelType OutputPixelType;

  typedef typename InputImageType::RegionType  InputImageRegionType;
  typedef typename OutputImageType::RegionType OutputImageRegionType;

  typedef typename InputImageType::SizeType InputSizeType;

  typedef itk::Image<unsigned char, 3>                            UCHARIMAGE;
  typedef itk::ImageMaskSpatialObject<UCHARIMAGE::ImageDimension> ImageMaskSpatialObjectType;

  /** */
  itkSetMacro(OtsuPercentileThreshold, double);
  itkGetConstMacro(OtsuPercentileThreshold, double);
  /** */
  itkSetMacro(ThresholdCorrectionFactor, double);
  itkGetConstMacro(ThresholdCorrectionFactor, double);
  /** The closing size in mm, this is rounded up to the next closest number of
    * voxel by taking Spacing into account */
  itkSetMacro(ClosingSize, double);
  itkGetConstMacro(ClosingSize, double);
  /** The dilation size in mm, this is rounded up to the next closest number of
    * voxel by taking Spacing into account */
  itkSetMacro(DilateSize, double);
  itkGetConstMacro(DilateSize, double);

  // NOTE:  This will generate a new spatial object each time it is called, and
  // not return the previous spatial object
  ImageMaskPointer GetSpatialObjectROI(void)
  {
    if( m_ResultMaskPointer.IsNull() )  // This is a cheap way to only create
                                        // the mask once, note that this is made
                                        // null when GenerateData is called.
      {
      typedef itk::CastImageFilter<OutputImageType, UCHARIMAGE> CastImageFilter;
      typename CastImageFilter::Pointer castFilter = CastImageFilter::New();
      castFilter->SetInput( this->GetOutput() );
      castFilter->Update();

      // convert mask image to mask
      typename ImageMaskSpatialObjectType::Pointer mask = ImageMaskSpatialObjectType::New();
      mask->SetImage( castFilter->GetOutput() );
      mask->ComputeObjectToWorldTransform();
      m_ResultMaskPointer = dynamic_cast<ImageMaskSpatialObjectType *>( mask.GetPointer() );
      if( m_ResultMaskPointer.IsNull() )
        {
        itkGenericExceptionMacro(<< "failed conversion to MaskSpatialObject");
        }
      }
    return m_ResultMaskPointer;
  }

  typename UCHARIMAGE::ConstPointer GetBinaryImageROI()
  {
    ImageMaskPointer tmp = this->GetSpatialObjectROI();

    typename UCHARIMAGE::ConstPointer rval = nullptr;
    if( tmp.IsNotNull() )
      {
      const typename itk::ImageMaskSpatialObject<3>::ConstPointer imso =
        dynamic_cast<itk::ImageMaskSpatialObject<3> *>(tmp.GetPointer() );
      if( imso.IsNull() )
        {
        itkGenericExceptionMacro(<< "failed conversion to MaskSpatialObject");
        }
      if( imso.IsNotNull() )
        {
        rval = imso->GetImage();
        }
      }
    return rval;
  }

#ifdef ITK_USE_CONCEPT_CHECKING
  /** Begin concept checking */
  itkConceptMacro( SameDimensionCheck,
                   ( Concept::SameDimension<InputImageDimension, OutputImageDimension> ) );
  /** End concept checking */
#endif
protected:
  BRAINSROIAutoImageFilter();
  virtual ~BRAINSROIAutoImageFilter()
  {
  }

  void PrintSelf(std::ostream & os, Indent indent) const ITK_OVERRIDE;

  void GenerateData() ITK_OVERRIDE;

private:
  ITK_DISALLOW_COPY_AND_ASSIGN(BRAINSROIAutoImageFilter);

  double           m_OtsuPercentileThreshold;
  double           m_ThresholdCorrectionFactor;
  double           m_ClosingSize;
  double           m_DilateSize;
  ImageMaskPointer m_ResultMaskPointer;
};
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkBRAINSROIAutoImageFilter.hxx"
#endif

#endif
