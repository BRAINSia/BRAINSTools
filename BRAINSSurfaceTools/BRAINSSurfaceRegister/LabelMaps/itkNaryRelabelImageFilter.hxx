/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkNaryRelabelImageFilter.txx,v $
  Language:  C++
  Date:      $Date: 2007/04/04 12:43:52 $
  Version:   $Revision: 1.21 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef _itkNaryRelabelImageFilter_txx
#define _itkNaryRelabelImageFilter_txx

#include "itkNaryRelabelImageFilter.h"
#include "itkImageRegionIterator.h"
#include "itkProgressReporter.h"

namespace itk
{
/**
 * Constructor
 */
template <class TInputImage, class TOutputImage>
NaryRelabelImageFilter<TInputImage, TOutputImage>
::NaryRelabelImageFilter()
{
  m_BackgroundValue = NumericTraits<InputImagePixelType>::Zero;
  m_IgnoreCollision = true;
}

/**
 * ThreadedGenerateData Performs the pixel-wise addition
 */
template <class TInputImage, class TOutputImage>
void
NaryRelabelImageFilter<TInputImage, TOutputImage>
::GenerateData()
{
  this->AllocateOutputs();

  const OutputImageRegionType & outputRegionForThread = this->GetOutput()->GetRequestedRegion();

  typedef ImageRegionConstIterator<TInputImage> ImageRegionConstIteratorType;
  std::vector<ImageRegionConstIteratorType *> inputIterators;
  // create the iterators for the input images
  for( unsigned int i = 0; i < this->GetNumberOfInputs(); ++i )
    {
    const InputImageType * input = this->GetInput( i );

    if( input )
      {
      ImageRegionConstIteratorType * inputIt = new ImageRegionConstIteratorType( input, outputRegionForThread );
      inputIterators.push_back( inputIt );
      }
    }

  // create the progress reporter
  ProgressReporter progress(this, 0, outputRegionForThread.GetNumberOfPixels() * 2 * inputIterators.size() );

  // found the labels in the input images and compute their new value
  typedef std::vector<std::map<InputImagePixelType, OutputImagePixelType> > TranslatorType;
  TranslatorType translator;
  translator.resize( inputIterators.size() );

  OutputImagePixelType label = NumericTraits<OutputImagePixelType>::Zero;
  for( unsigned int i = 0; i < inputIterators.size(); ++i )
    {
    ImageRegionConstIteratorType * inputIt = inputIterators[i];
    for( inputIt->GoToBegin(); !inputIt->IsAtEnd(); ++(*inputIt) )
      {
      const InputImagePixelType & v = inputIt->Get();
      if( v != m_BackgroundValue && translator[i].find( v ) == translator[i].end() )
        {
        // a new label to translate
        if( label == m_BackgroundValue )
          {
          // avoid the background label
          label++;
          }
        // std::cout << i << " " << v+0.0 << " " << label+0.0 << std::endl;
        translator[i][v] = label;

        // increment the label for the next to translate
        // TODO: throw an exception if the maximum number of labels is exceeded
        label++;
        }
      progress.CompletedPixel();
      }
    }

  // now write the output image
  OutputImagePointer                output = this->GetOutput();
  ImageRegionIterator<TOutputImage> outputIt( output, outputRegionForThread );
  for( unsigned int i = 0; i < inputIterators.size(); ++i )
    {
    inputIterators[i]->GoToBegin();
    }
  for( outputIt.GoToBegin(); !outputIt.IsAtEnd(); ++outputIt )
    {
    // find the output label
    OutputImagePixelType outputLabel = m_BackgroundValue;
    bool                 labelFound = false;
    for( unsigned int i = 0; i < inputIterators.size(); i++ )
      {
      const InputImagePixelType & v = inputIterators[i]->Get();
      if( v != m_BackgroundValue )
        {
        if( !m_IgnoreCollision && labelFound )
          {
          itkExceptionMacro( << "Label collision detected." );
          }

        outputLabel = translator[i][v];
        labelFound = true;
        }
      }

    // write the output
    outputIt.Set( outputLabel );
    // now update the input iterators
    for( unsigned int i = 0; i < inputIterators.size(); i++ )
      {
      ++(*inputIterators[i]);
      progress.CompletedPixel();
      }
    }
  // delete the input iterators
  for( unsigned int i = 0; i < inputIterators.size(); ++i )
    {
    delete inputIterators[i];
    }
}

template <class TInputImage, class TOutputImage>
void
NaryRelabelImageFilter<TInputImage, TOutputImage>
::PrintSelf( std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf( os, indent );

  os << indent << "Background Value: "
     << static_cast<typename NumericTraits<InputImagePixelType>::PrintType>(m_BackgroundValue) << std::endl;
  os << indent << "Ignore Collision: " << static_cast<typename NumericTraits<bool>::PrintType>(m_IgnoreCollision)
     << std::endl;
}
} // end namespace itk

#endif
