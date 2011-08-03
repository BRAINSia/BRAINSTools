/*
 *  itkResampleInPlaceImageFilter.hxx
 *
 *
 *  Created by Wei Lu on 10/14/10.
 *
 */

#ifndef __itkResampleInPlaceImageFilter_hxx
#define __itkResampleInPlaceImageFilter_hxx

#include "itkResampleInPlaceImageFilter.h"
#include "itkCastImageFilter.h"

namespace itk
{
/**
 * Constructor
 */
template <class TInputImage, class TOutputImage>
ResampleInPlaceImageFilter<TInputImage, TOutputImage>
::ResampleInPlaceImageFilter() :
  m_OutputImage( NULL ),
  m_RigidTransform( NULL )
{
  this->SetNumberOfRequiredInputs( 1 );
}

/**
 * Set/Get input image, required
 */
template <class TInputImage, class TOutputImage>
void
ResampleInPlaceImageFilter<TInputImage, TOutputImage>
::SetInputImage( const InputImageType * image )
{
  this->SetInput( 0, image );
}

template <class TInputImage, class TOutputImage>
const typename ResampleInPlaceImageFilter<TInputImage, TOutputImage>::InputImageType
* ResampleInPlaceImageFilter<TInputImage, TOutputImage>
::GetInputImage() const
  {
  return this->GetInput( 0 );
  }

/**
 * GenerateData Performs the in-place resampling
 */
template <class TInputImage, class TOutputImage>
void
ResampleInPlaceImageFilter<TInputImage, TOutputImage>
::GenerateData()
{
  if( !this->GetInput() )
    {
    itkExceptionMacro(<< "Input image has not been connected");
    return;
    }

  typedef typename RigidTransformType::Pointer RigidTransformPointer;
  RigidTransformPointer invOfRigidTransform = RigidTransformType::New();
  invOfRigidTransform->SetIdentity();
  const InputImagePointType centerPoint = m_RigidTransform->GetCenter();
  invOfRigidTransform->SetCenter( centerPoint );
  invOfRigidTransform->SetIdentity();
  this->m_RigidTransform->GetInverse( invOfRigidTransform );    // Cache the
                                                                // inverse

    {
    /** make a cast copied version of the input image **/
    typedef CastImageFilter<InputImageType, OutputImageType> DuplicatorType;
    typename DuplicatorType::Pointer CastFilter = DuplicatorType::New();
    CastFilter->SetInput( this->GetInput() );
    CastFilter->Update();
    m_OutputImage = CastFilter->GetOutput();
    }

  // Modify the origin and direction info of the image to reflect the transform.
  m_OutputImage->SetOrigin( invOfRigidTransform->GetMatrix()
                            * this->GetInput()->GetOrigin() + invOfRigidTransform->GetTranslation() );
  m_OutputImage->SetDirection( invOfRigidTransform->GetMatrix()
                               * this->GetInput()->GetDirection() );

  this->GraftOutput( m_OutputImage );
}

template <class TInputImage, class TOutputImage>
void
ResampleInPlaceImageFilter<TInputImage, TOutputImage>::PrintSelf( std::ostream& os, Indent indent ) const
{
  Superclass::PrintSelf( os, indent );

  os << indent << "Input transform: " << m_RigidTransform << std::endl;
  os << indent << "Output image: " << m_OutputImage << std::endl;
}
} // end namespace itk

#endif
