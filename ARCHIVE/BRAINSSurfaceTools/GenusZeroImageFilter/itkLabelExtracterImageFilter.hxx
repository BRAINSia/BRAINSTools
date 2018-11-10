/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkLabelExtracterImageFilter.hxx,v $
  Language:  C++
  Date:      $Date: 2009/07/17 14:37:03 $
  Version:   $Revision: 1.1 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

  Portions of this code are covered under the VTK copyright.
  See VTKCopyright.txt or http://www.kitware.com/VTKCopyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef _itkLabelExtracterImageFilter_hxx
#define _itkLabelExtracterImageFilter_hxx

#include "itkLabelExtracterImageFilter.h"

namespace itk
{
/**
 *
 */
template <typename TInputImage, typename TOutputImage>
LabelExtracterImageFilter<TInputImage, TOutputImage>
::LabelExtracterImageFilter()
{
}

/**
 *
 */
template <typename TInputImage, typename TOutputImage>
void
LabelExtracterImageFilter<TInputImage, TOutputImage>
::SetChange( const InputPixelType & original, const OutputPixelType & result )
{
  OutputPixelType current = this->GetFunctor().GetChange(original);

  if( current != result )
    {
    this->GetFunctor().SetChange(original, result);
    this->Modified();
    }
}

/**
 *
 */
template <typename TInputImage, typename TOutputImage>
void
LabelExtracterImageFilter<TInputImage, TOutputImage>
::SetChangeMap( const ChangeMapType & changeMap )
{
  // If the whole map is being set then we assume that a real change is made
  this->GetFunctor().SetChangeMap(changeMap);
  this->Modified();
}

/**
 *
 */
template <typename TInputImage, typename TOutputImage>
void
LabelExtracterImageFilter<TInputImage, TOutputImage>
::ClearChangeMap()
{
  // If the whole map is being set then we assume that a real change is made
  this->GetFunctor().ClearChangeMap();
  this->Modified();
}

/**
 *
 */
template <typename TInputImage, typename TOutputImage>
void
LabelExtracterImageFilter<TInputImage, TOutputImage>
::PrintSelf(std::ostream & os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);
  // Maybe should iterate the change map and print it here
}

} // end namespace itk

#endif
