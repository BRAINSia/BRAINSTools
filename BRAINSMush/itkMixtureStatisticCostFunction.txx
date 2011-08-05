/*=========================================================================
 *
 *  Program:   Insight Segmentation & Registration Toolkit
 *  Module:    $RCSfile: itkMixtureStatisticCostFunction.txx,v $
 *  Language:  C++
 *  Date:      $Date: 2007-03-29 19:37:00 $
 *  Version:   $Revision: 1.14 $
 *
 *  Copyright (c) Insight Software Consortium. All rights reserved.
 *  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.
 *
 *  This software is distributed WITHOUT ANY WARRANTY; without even
 *  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 *  PURPOSE.  See the above copyright notices for more information.
 *
 *  =========================================================================*/
#ifndef _itkMixtureStatisticCostFunction_txx
#define _itkMixtureStatisticCostFunction_txx

#include "itkMixtureStatisticCostFunction.h"
#include "assert.h"

namespace itk
{
template <class TFirstImage, class TSecondImage>
MixtureStatisticCostFunction<TFirstImage, TSecondImage>
::MixtureStatisticCostFunction()
{
  m_MeasurePointer = new MeasureType();
}

template <class TFirstImage, class TSecondImage>
MixtureStatisticCostFunction<TFirstImage, TSecondImage>
::~MixtureStatisticCostFunction()
{
  delete m_MeasurePointer;
}

template <class TFirstImage, class TSecondImage>
typename MixtureStatisticCostFunction<TFirstImage, TSecondImage>::MeasureType
MixtureStatisticCostFunction<TFirstImage, TSecondImage>
::GetValue(const ParametersType & parameters) const
{
  const double firstImageWeighting = parameters[0];
  const double secondImageWeighting = parameters[1];

  double num = m_NumberOfMaskVoxels;
  // weighted sum of additive statistics
  const double firstsum = firstImageWeighting * m_SumOfFirstMaskVoxels;
  const double secondsum = secondImageWeighting * m_SumOfSecondMaskVoxels;
  const double firstsumsq = firstImageWeighting * firstImageWeighting
    * m_SumSquaresOfFirstMaskVoxels;
  const double secondsumsq = secondImageWeighting * secondImageWeighting
    * m_SumSquaresOfSecondMaskVoxels;
  const double crosstermsum = 2.0 * firstImageWeighting
    * secondImageWeighting
    * m_SumOfFirstTimesSecondMaskVoxels;
  double sum = firstsum + secondsum;
  double sumsq = firstsumsq + secondsumsq + crosstermsum;

  // convert 3 statistics into mixture mean and variance
  sumsq = ( sumsq - ( sum * sum ) / num ) / ( num - 1 );
  sum = sum / num;
  // the measure is the squared distance to each of the goals.
  sum = sum - m_DesiredMean;
  sumsq = sumsq - m_DesiredVariance;
  m_Measure[0] = sum * sum;
  m_Measure[1] = sumsq * sumsq;
  return m_Measure;
}

template <class TFirstImage, class TSecondImage>
typename MixtureStatisticCostFunction<TFirstImage, TSecondImage>::MeasureType
* MixtureStatisticCostFunction<TFirstImage, TSecondImage>
::GetValue(ParametersType & parameters)
  {
  ( *m_MeasurePointer ) = GetValue(parameters);
  return m_MeasurePointer;
  }

template <class TFirstImage, class TSecondImage>
unsigned int
MixtureStatisticCostFunction<TFirstImage, TSecondImage>
::GetNumberOfParameters() const
{
  // Return the number of parameters.
  return 2;
}

template <class TFirstImage, class TSecondImage>
unsigned int
MixtureStatisticCostFunction<TFirstImage, TSecondImage>
::GetNumberOfValues() const
{
  // Return the number of residuals.
  return 2;
}

template <class TFirstImage, class TSecondImage>
void
MixtureStatisticCostFunction<TFirstImage, TSecondImage>
::Initialize(short label)
{
  m_Measure.SetSize(2);
  m_MeasurePointer->SetSize(2);

  // measure each image and each squared image within the mask

  typedef typename itk::ImageRegionConstIterator<typename Self::FirstImageType>
    FirstConstIteratorType;
  FirstConstIteratorType firstIt( m_FirstImage,
                                  m_FirstImage->GetRequestedRegion() );

  typedef typename itk::ImageRegionConstIterator<typename Self::SecondImageType>
    SecondConstIteratorType;
  SecondConstIteratorType secondIt( m_SecondImage,
                                    m_SecondImage->GetRequestedRegion() );

  typedef typename itk::ImageRegionConstIterator<typename Self::ImageMaskType>
    MaskConstIteratorType;
  MaskConstIteratorType maskIt( m_ImageMask, m_ImageMask->GetRequestedRegion() );

  m_NumberOfMaskVoxels = 0.0;
  m_SumOfFirstMaskVoxels = 0.0;
  m_SumOfSecondMaskVoxels = 0.0;
  m_SumSquaresOfFirstMaskVoxels = 0.0;
  m_SumSquaresOfSecondMaskVoxels = 0.0;
  m_SumOfFirstTimesSecondMaskVoxels = 0.0;
  for( maskIt.GoToBegin(), firstIt.GoToBegin(), secondIt.GoToBegin();
       !maskIt.IsAtEnd();
       ++maskIt, ++firstIt, ++secondIt )
    {
    if( maskIt.Get() == label )
      {
      FirstImagePixelType  firstValue = firstIt.Get();
      SecondImagePixelType secondValue = secondIt.Get();

      m_NumberOfMaskVoxels += 1.0;
      m_SumOfFirstMaskVoxels += firstValue;
      m_SumOfSecondMaskVoxels += secondValue;
      m_SumSquaresOfFirstMaskVoxels += firstValue * firstValue;
      m_SumSquaresOfSecondMaskVoxels += secondValue * secondValue;
      m_SumOfFirstTimesSecondMaskVoxels += firstValue * secondValue;
      }
    }
}

template <class TFirstImage, class TSecondImage>
void
MixtureStatisticCostFunction<TFirstImage, TSecondImage>
::PrintSelf(std::ostream & os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);

  os << indent << "DesiredMean = " << m_DesiredMean << std::endl
     << indent << "DesiredVariance = " << m_DesiredVariance << std::endl
     << indent << "NumberOfMaskVoxels = " << m_NumberOfMaskVoxels << std::endl
     << indent << "SumOfFirstMaskVoxels = " << m_SumOfFirstMaskVoxels
     << std::endl
     << indent << "SumSquaresOfFirstMaskVoxels = "
     << m_SumSquaresOfFirstMaskVoxels << std::endl
     << indent << "SumOfSecondMaskVoxels = " << m_SumOfSecondMaskVoxels
     << std::endl
     << indent << "SumSquaresOfSecondMaskVoxels = "
     << m_SumSquaresOfSecondMaskVoxels << std::endl;
}
}   // end namespace itk
#endif
