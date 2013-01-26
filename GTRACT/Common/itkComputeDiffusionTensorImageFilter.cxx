/*=========================================================================

 Program:   GTRACT (Guided Tensor Restore Anatomical Connectivity Tractography)
 Module:    $RCSfile: $
 Language:  C++
 Date:      $Date: 2006/03/29 14:53:40 $
 Version:   $Revision: 1.9 $

   Copyright (c) University of Iowa Department of Radiology. All rights reserved.
   See GTRACT-Copyright.txt or http://mri.radiology.uiowa.edu/copyright/GTRACT-Copyright.txt
   for details.

      This software is distributed WITHOUT ANY WARRANTY; without even
      the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
      PURPOSE.  See the above copyright notices for more information.

=========================================================================*/

#ifndef __itkComputeDiffusionTensorImageFilter_cxx
#define __itkComputeDiffusionTensorImageFilter_cxx

#include "itkImageRegionIterator.h"
#include "itkImageRegionConstIterator.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkImageLinearConstIteratorWithIndex.h"
#include "itkConstNeighborhoodIterator.h"
#include "itkNeighborhoodAlgorithm.h"
#include <itkIOCommon.h>
#include "itkMetaDataObject.h"
#include "itkProgressAccumulator.h"
#include "itkMedianImageFilter.h"

#include "itkComputeDiffusionTensorImageFilter.h"
#include "algo.h"

#include <iostream>

namespace itk
{
ComputeDiffusionTensorImageFilter
::ComputeDiffusionTensorImageFilter()
{
  m_UseMedianFilter     = false;
  m_BackgroundThreshold = 0;
  m_NumberOfDirections  = 0;
  m_NumberOfBSteps      = 0;
  m_MedianFilterSize.Fill(0);
}

void
ComputeDiffusionTensorImageFilter
::Update()
{
  if( m_UseMedianFilter == true )
    {
    typedef itk::MedianImageFilter<InputImageType, InputImageType> MedianFilterType;
    MedianFilterType::Pointer filter = MedianFilterType::New();
    filter->SetInput( m_Input );
    filter->SetRadius( m_MedianFilterSize );
    filter->Update();
    m_InternalImage = filter->GetOutput();
    }
  else
    {
    m_InternalImage = m_Input;
    }

  std::cout << "Tensor Directions: " << std::endl;
  std::cout << m_DiffusionDirections << std::endl;

  OutputImageRegionType  TensorRegion;
  OutputImageSizeType    TensorSize;
  OutputImageIndexType   TensorIndex;
  OutputImageSpacingType TensorSpacing;
  OutputImagePointType   TensorOrigin;

  InputImageRegionType  ADCRegion  = m_InternalImage->GetLargestPossibleRegion();
  InputImageSizeType    ADCSize    = ADCRegion.GetSize();
  InputImageIndexType   ADCIndex  = ADCRegion.GetIndex();
  InputImageSpacingType ADCSpacing  = m_InternalImage->GetSpacing();
  InputImagePointType   ADCOrigin  = m_InternalImage->GetOrigin();
  for( int i = 0; i < 3; i++ )
    {
    TensorSize[i]     = ADCSize[i];
    TensorIndex[i]   = ADCIndex[i];
    TensorSpacing[i] = ADCSpacing[i];
    TensorOrigin[i]  = ADCOrigin[i];
    }

  TensorRegion.SetSize( TensorSize );
  TensorRegion.SetIndex( TensorIndex );

  m_Output = OutputImageType::New();
  m_Output->SetRegions( TensorRegion );
  m_Output->SetSpacing( TensorSpacing );
  m_Output->SetOrigin( TensorOrigin );
  m_Output->Allocate();

  typedef itk::ImageLinearConstIteratorWithIndex<InputImageType> ADCIteratorType;
  ADCIteratorType ADC_It( m_InternalImage, m_InternalImage->GetRequestedRegion() );

  ADC_It.SetDirection(3);
  ADC_It.GoToBegin();
  TMatrix mMatrix = Matrix_Inverse(m_DiffusionDirections);

  while( !ADC_It.IsAtEnd() )
    {
    /* Get the index & value for B0 */
    ADC_It.GoToBeginOfLine();
    float ADC0 = (float)ADC_It.Get();
    ADCIndex = ADC_It.GetIndex();
    for( int i = 0; i < 3; i++ )
      {
      TensorIndex[i] = ADCIndex[i];
      }

    OutputPixelType currentVoxel;
    currentVoxel.Fill(0);

    if( ADC0 > m_BackgroundThreshold )
      {
      ++ADC_It;
      TVector temp(m_NumberOfBSteps * m_NumberOfDirections);
      for( int i = 0; i < m_NumberOfBSteps * m_NumberOfDirections; i++ )
        {
        temp(i) = (float)ADC_It.Get();
        ++ADC_It;
        }

      // ////////////////////////////////////////////////////////////////////////
      TVector ADCm(m_NumberOfDirections);
      TVector Ln_ADCs(m_NumberOfBSteps + 1);
      Ln_ADCs(0) = 0;
      bool ErrFlg = false;
      for( int direction = 0; direction < m_NumberOfDirections; direction++ )
        {
        for( int step = 0; step < m_NumberOfBSteps; step++ )
          {
          float tempflt;
          tempflt = temp(direction * m_NumberOfBSteps + step);

          if( tempflt == 0 )
            {
            ErrFlg = true;  break;
            }
          Ln_ADCs(step + 1) = vcl_log(tempflt / ADC0);
          }
        if( ErrFlg )
          {
          break;
          }
        ADCm(direction) = -1 * My_lsf(m_BValues, Ln_ADCs);
        }
      if( !ErrFlg )
        {
        TVector ADCe = mMatrix * ADCm;
        currentVoxel.SetVnlVector(ADCe);
        }
      m_Output->SetPixel(TensorIndex, currentVoxel);
      }
    else
      {
      m_Output->SetPixel(TensorIndex, currentVoxel);
      }
    // ////////////////////////////////////////////////////////////////////////

    ADC_It.NextLine();
    }

  itk::Point<double, 3> fixedOrigin = m_Output->GetOrigin();
  fixedOrigin.GetVnlVector().fill(0.0);
  m_Output->SetOrigin(fixedOrigin);
  m_Output->SetMetaDataDictionary( m_Input->GetMetaDataDictionary() );
}
} // end namespace itk
#endif
