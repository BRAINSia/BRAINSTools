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

#ifndef __itkFixDicomImageFilter_hxx
#define __itkFixDicomImageFilter_hxx

#include "itkImageRegionIterator.h"
#include "itkImageRegionConstIterator.h"
#include "itkOrientImageFilter.h"
#include "itkFixDicomImageFilter.h"
#include <itkIOCommon.h>
#include <itkCastImageFilter.h>
#include <itkConstantPadImageFilter.h>
#include <itkExtractImageFilter.h>
#include "itkMetaDataObject.h"
#include "itkProgressAccumulator.h"

#include <iostream>

namespace itk
{
template <class TInputImage, class TOutputImage>
FixDicomImageFilter<TInputImage, TOutputImage>
::FixDicomImageFilter()
{
  m_NumberOfImagesPerSlice = 1;
  m_NumberBValues = 0;
  m_NumberDtiDirections = 0;
  m_DirectionsIncrementFirst = true;
}

template <class TInputImage, class TOutputImage>
void
FixDicomImageFilter<TInputImage, TOutputImage>
::SetDirectionsIncrementFirst()
{
  m_DirectionsIncrementFirst = true;
}

template <class TInputImage, class TOutputImage>
void
FixDicomImageFilter<TInputImage, TOutputImage>
::SetBvaluesIncrementFirst()
{
  m_DirectionsIncrementFirst = false;
}

template <class TInputImage, class TOutputImage>
void
FixDicomImageFilter<TInputImage, TOutputImage>
::Update()
{
  InputImageRegionType    region3D  = m_Input->GetLargestPossibleRegion();
  InputImageIndexType     index3D    = region3D.GetIndex();
  InputImageSizeType      size3D    = region3D.GetSize();
  InputImageSpacingType   spacing3D  = m_Input->GetSpacing();
  InputImagePointType     origin3D  = m_Input->GetOrigin();
  InputImageDirectionType direction3D = m_Input->GetDirection();

  OutputImageRegionType    region4D;
  OutputImageIndexType     index4D;
  OutputImageSizeType      size4D;
  OutputImageSpacingType   spacing4D;
  OutputImagePointType     origin4D;
  OutputImageDirectionType direction4D;

  for( int i = 0; i < 3; i++ )
    {
    index4D[i]  = index3D[i];
    size4D[i]  = size3D[i];
    spacing4D[i] = spacing3D[i];
    origin4D[i] = origin3D[i];
    for( int j = 0; j < 3; j++ )
      {
      direction4D[i][j] = direction3D[i][j];
      }
    }
  index4D[3]   = 0;
  spacing4D[3] = 1;
  origin4D[3] = 0;
  size4D[3] = m_NumberOfImagesPerSlice;
  direction4D[3][0] = 0.0; direction4D[3][1] = 0.0; direction4D[3][2] = 0.0; direction4D[3][3] = 1.0;

  const int numberOfSlice = int(size3D[2] / size4D[3]);
  size4D[2] = numberOfSlice;
  region4D.SetIndex(index4D);
  region4D.SetSize(size4D);

  m_Output = OutputImageType::New();
  m_Output->SetRegions(region4D);
  m_Output->SetSpacing(spacing4D);
  m_Output->SetOrigin(origin4D);
  m_Output->Allocate();

  /* Set MataData Dictionary and Direction Cosines*/
  m_Output->SetMetaDataDictionary( m_Input->GetMetaDataDictionary() );
  m_Output->SetDirection(direction4D);

  typedef itk::ImageLinearConstIteratorWithIndex<TInputImage> IteratorType;
  IteratorType it( m_Input, m_Input->GetRequestedRegion() );

  it.SetDirection(2);
  it.GoToBegin();
  while( !it.IsAtEnd() )
    {
    it.GoToBeginOfLine();
    while( !it.IsAtEndOfLine() )
      {
      index3D = it.GetIndex();
      index4D[0] = index3D[0]; index4D[1] = index3D[1];
      index4D[2] = index3D[2] % numberOfSlice;
      if( m_NumberBValues > 0 )
        {
        int tmpIndex = int(index3D[2] / numberOfSlice);
        if( tmpIndex == 0 )
          {
          index4D[3] = 0;
          }
        else
          {
          if( m_DirectionsIncrementFirst )
            {
            int tmpOffset = ( tmpIndex - 1 ) % m_NumberDtiDirections;
            int tmpScale =  ( tmpIndex - 1 ) / m_NumberDtiDirections;
            index4D[3] = tmpOffset * m_NumberBValues + tmpScale + 1;
            }
          else
            {
            index4D[3] = tmpIndex;
            }
          }
        }
      else
        {
        index4D[3] = int(index3D[2] / numberOfSlice);
        }
      m_Output->SetPixel( index4D, it.Get() );
      ++it;
      }

    it.NextLine();
    }
}
} // end namespace itk
#endif
