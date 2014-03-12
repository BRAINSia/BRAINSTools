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

#ifndef __itkOrient4dImageFilter_hxx
#define __itkOrient4dImageFilter_hxx

#include "itkImageRegionIterator.h"
#include "itkImageRegionConstIterator.h"
#include "itkOrientImageFilter.h"
#include "itkOrient4dImageFilter.h"
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
Orient4dImageFilter<TInputImage, TOutputImage>
::Orient4dImageFilter()
{
  m_ExtractImageFilter = ExtractFilterType::New();
  m_OrientImageFilter  = OrientFilterType::New();

  m_FlipXaxis = 0;
  m_FlipYaxis = 0;
  m_FlipZaxis = 0;
}

template <class TInputImage, class TOutputImage>
void
Orient4dImageFilter<TInputImage, TOutputImage>
::Update()
{
  InputImageRegionType  imageRegion  = m_Input->GetLargestPossibleRegion();
  InputImageSizeType    imageSize    = imageRegion.GetSize();
  InputImageSpacingType imageSpacing = m_Input->GetSpacing();
  InputImagePointType   imageOrigin  = m_Input->GetOrigin();
  const int             numVolumes = imageSize[3];

  // std::cout << "Volumes : " << numVolumes << std::endl;

  imageSize[3] = 0;
  imageRegion.SetSize(imageSize);

  m_Output = TOutputImage::New();
  /*
  m_Output->SetRegions(imageRegion);
  m_Output->SetSpacing(imageSpacing);
  m_Output->SetOrigin(imageOrigin);
  m_Output->Allocate();
   */
  OutputImageIndexType imageIndex = imageRegion.GetIndex();
  imageIndex[0] = 0;
  imageIndex[1] = 0;
  imageIndex[2] = 0;

  // std::cout << "Input Image : " << m_Input << std::endl;
  int /*xsize,*/ ysize, zsize;
  for( int i = 0; i < numVolumes; i++ )
    {
    // std::cout << "Orient : " << i << std::endl;
    imageIndex[3] = i;
    imageRegion.SetIndex(imageIndex);
    // std::cout << "Extract : " << imageRegion << std::endl;
    m_ExtractImageFilter->SetExtractionRegion( imageRegion );
    m_ExtractImageFilter->SetInput(m_Input);
    m_ExtractImageFilter->Update();
    // std::cout << "Extract Region : " << imageIndex << std::endl;

    ExtractImagePointer Extract3DImage = m_ExtractImageFilter->GetOutput();

    // std::cout << "Set Dirs : " << std::endl;

    if( m_OrientImageFilter->GetUseImageDirection() )
      {
      InputImageDirectionType   directions4d = m_Input->GetDirection();
      ExtractImageDirectionType directions3d = Extract3DImage->GetDirection();
      // std::cout << "Original Directions : " << directions4d << std::endl;
      for( int n = 0; n < 3; n++ )
        {
        for( int m = 0; m < 3; m++ )
          {
          directions3d[n][m] = directions4d[n][m];
          }
        }
      Extract3DImage->SetDirection( directions3d );
      // std::cout << "3d Directions : " << directions3d << std::endl;
      }

    // std::cout << "Orient : " << std::endl;
    m_OrientImageFilter->SetInput( Extract3DImage );

    m_OrientImageFilter->Update();

    FlipFilterPointerType flipImageFilter = FlipFilterType::New();
    FlipFilterAxesType    flipAxis = flipImageFilter->GetFlipAxes();
    flipAxis[0] = m_FlipXaxis;
    flipAxis[1] = m_FlipYaxis;
    flipAxis[2] = m_FlipZaxis;
    flipImageFilter->SetInput( m_OrientImageFilter->GetOutput() );
    flipImageFilter->SetFlipAxes( flipAxis );
    flipImageFilter->Update();

    OrientImagePointer Orient3DImage = flipImageFilter->GetOutput();

    // std::cout << "Copy Image : " << std::endl;
    if( i == 0 )
      {
      InputImageRegionType  newRegion  = m_Output->GetLargestPossibleRegion();
      InputImageSpacingType newSpacing = m_Output->GetSpacing();
      InputImagePointType   newOrigin  = m_Output->GetOrigin();
      InputImageSizeType    newSize    = newRegion.GetSize();

      OrientImageRegionType  orRegion  = Orient3DImage->GetLargestPossibleRegion();
      OrientImageSpacingType orSpacing = Orient3DImage->GetSpacing();
      OrientImagePointType   orOrigin  = Orient3DImage->GetOrigin();
      OrientImageSizeType    orSize    = orRegion.GetSize();
      for( int m = 0; m < 3; m++ )
        {
        newSpacing[m] = orSpacing[m];
        newOrigin[m] = orOrigin[m];
        newSize[m] = orSize[m];
        }
      // xsize =   newSize[0];
      ysize =   newSize[1];
      zsize =   newSize[2];
      newSize[3] = numVolumes;
      newRegion.SetSize(newSize);
      newSpacing[3] = imageSpacing[3];
      newOrigin[3] = imageOrigin[3];
      // std::cout << "Init Image : " << std::endl;
      m_Output->SetRegions(newRegion);
      m_Output->SetSpacing(newSpacing);
      m_Output->SetOrigin(newOrigin);
      m_Output->Allocate();
      }
    for( int z = 0; z < zsize; z++ )
      {
      for( int y = 0; y < ysize; y++ )
        {
        for( int x = 0; x < zsize; x++ )
          {
          OutputImageIndexType index4D;
          OrientImageIndexType index3D;
          index3D[0] = index4D[0] = x;
          index3D[1] = index4D[1] = y;
          index3D[2] = index4D[2] = z;
          index4D[3] = i;
          m_Output->SetPixel( index4D, Orient3DImage->GetPixel(index3D) );
          }
        }
      }
    }

  // Set Meta Data Orientation Information

  itk::MetaDataDictionary imageDictionary;
  itk::EncapsulateMetaData<itk::SpatialOrientation::ValidCoordinateOrientationFlags>(
    imageDictionary,
    ITK_CoordinateOrientation,
    m_OrientImageFilter->
    GetDesiredCoordinateOrientation() );

  m_Output->SetMetaDataDictionary(imageDictionary);
}
} // end namespace itk
#endif
