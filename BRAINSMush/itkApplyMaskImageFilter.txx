#ifndef _itkApplyMaskImageFilter_txx
#define _itkApplyMaskImageFilter_txx

#include "itkApplyMaskImageFilter.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionConstIterator.h"

namespace itk
{
template <class TInputImage, class TOutputImage>
ApplyMaskImageFilter<TInputImage, TOutputImage>
::ApplyMaskImageFilter()
{
  this->SetNumberOfRequiredInputs(2);

  m_InvertMask = false;
}

template <class TInputImage, class TOutputImage>
void
ApplyMaskImageFilter<TInputImage, TOutputImage>
::PrintSelf(std::ostream & os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);

  os << indent << "InvertMask: ";
  os << m_InvertMask << std::endl;
}

template <class TInputImage, class TOutputImage>
void
ApplyMaskImageFilter<TInputImage, TOutputImage>
::SetMaskImage(const InputImageType *reference)
{
  this->ProcessObject::SetNthInput( 1,
                                    const_cast<InputImageType *>( reference ) );
}

template <class TInputImage, class TOutputImage>
const typename ApplyMaskImageFilter<TInputImage, TOutputImage>
::InputImageType
* ApplyMaskImageFilter<TInputImage, TOutputImage>
::GetMaskImage()
  {
  if( this->GetNumberOfInputs() < 2 )
    {
    return NULL;
    }

  return dynamic_cast<TInputImage *>(
    this->ProcessObject::GetInput(1) );
  }

template <class TInputImage, class TOutputImage>
void
ApplyMaskImageFilter<TInputImage, TOutputImage>
::GenerateData()
{
  this->AllocateOutputs();

  // get the input and output pointers
  InputImageConstPointer inputImage  = this->GetInput(0);
  InputImageConstPointer maskImage   = this->GetInput(1);
  OutputImagePointer     outputVolume = this->GetOutput();

  // mask the image
  typedef ImageRegionConstIterator<InputImageType> InputConstIterator;
  typedef ImageRegionIterator<OutputImageType>     OutputIterator;

  InputConstIterator itImage( inputImage,
                              inputImage->GetRequestedRegion() );

  InputConstIterator itMask( maskImage,
                             maskImage->GetRequestedRegion() );

  OutputIterator itOut( outputVolume,
                        outputVolume->GetRequestedRegion() );

  // support progress methods/callbacks
  unsigned long updateVisits = 0;
  unsigned long totalPixels  = 0;

  totalPixels  = inputImage->GetRequestedRegion().GetNumberOfPixels();
  updateVisits = totalPixels / 10;
  if( updateVisits < 1 )
    {
    updateVisits = 1;
    }

  itImage.GoToBegin();
  itMask.GoToBegin();
  unsigned long i = 0;
  for( itOut.GoToBegin(); !itOut.IsAtEnd(); ++itOut )
    {
    if( !( i % updateVisits ) )
      {
      this->UpdateProgress( (float)i / (float)totalPixels );
      }

    if( !m_InvertMask )
      {
      if( itMask.Get() != 0 )
        {
        itOut.Set( static_cast<OutputPixelType>( itImage.Get() ) );
        }
      else
        {
        itOut.Set( static_cast<OutputPixelType>( 0 ) );
        }
      }
    else
      {
      if( itMask.Get() != 0 )
        {
        itOut.Set( static_cast<OutputPixelType>( 0 ) );
        }
      else
        {
        itOut.Set( static_cast<OutputPixelType>( itImage.Get() ) );
        }
      }
    i++;
    ++itMask;
    ++itImage;
    }
}
}   // end namespace itk

#endif
