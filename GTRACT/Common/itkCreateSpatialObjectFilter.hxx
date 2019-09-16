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

#ifndef __itkCreateSpatialObjectFilter_hxx
#define __itkCreateSpatialObjectFilter_hxx

#include "itkImageRegionIterator.h"
#include "itkImageRegionConstIterator.h"
#include "itkOrientImageFilter.h"
#include "itkCreateSpatialObjectFilter.h"
#include <itkIOCommon.h>
#include <itkCastImageFilter.h>
#include <itkConstantPadImageFilter.h>
#include <itkExtractImageFilter.h>
#include "itkMetaDataObject.h"
#include "itkProgressAccumulator.h"

#include <iostream>

namespace itk
{
template <typename TInputImage, typename TTransformType, typename TSpatialObject>
CreateSpatialObjectFilter<TInputImage, TTransformType, TSpatialObject>::CreateSpatialObjectFilter()
{
  m_Input = NULL;
  m_Output = NULL;
  m_Transform = NULL;
}

template <typename TInputImage, typename TTransformType, typename TSpatialObject>
void
CreateSpatialObjectFilter<TInputImage, TTransformType, TSpatialObject>::Update()
{
  SpatialObjectPointListType points;

  m_Output = SpatialObjectType::New();

  using IteratorType = itk::ImageRegionConstIteratorWithIndex<InputImageType>;
  IteratorType it(m_Input, m_Input->GetLargestPossibleRegion());
  // int count;
  for (it.GoToBegin(); !it.IsAtEnd(); ++it)
  {
    if (it.Get() > 0)
    {
      //      count++;
      InputImageIndexType index = it.GetIndex();
      InputImagePointType p;
      m_Input->TransformIndexToPhysicalPoint(index, p);
      p = m_Transform->TransformPoint(p);
      BlobPointType blobPoint;
      blobPoint.SetPosition(p);
      points.push_back(blobPoint);
    }
  }

  m_Output->SetPoints(points);
  m_Output->ComputeBoundingBox();

  return;
}

/*
InputImagePointer maskImage;

if (m_Transform.IsNotNull())
{
using ResampleFilterType = itk::ResampleImageFilter<
InputImageType,
InputImageType >;

ResampleFilterType::Pointer resampler = ResampleFilterType::New();
//TransformPointer inverseTransform = TTransformType::New();
//inverseTransform->SetCenter(m_Transform->GetCenter());
//m_Transform->GetInverse( inverseTransform );
resampler->SetTransform( m_Transform );
resampler->SetInput( m_Image );
resampler->SetSize(    m_Size );
resampler->SetOutputOrigin(  m_Origin );
resampler->SetOutputSpacing( m_Spacing );
resampler->SetOutputDirection( m_Direction );
resampler->SetDefaultPixelValue( 0 ); // DEBUG
resampler->Update();

maskImage = resampler->GetOutput();
}
else
{
maskImage = m_Image;
}



using IteratorType = itk::ImageRegionConstIteratorWithIndex< InputImageType >;
IteratorType it(maskImage, maskImage->GetLargestPossibleRegion());

using BlobPointType = itk::SpatialObjectPoint<3>;
SpatialObjectType::PointListType    points;
m_Blob  = SpatialObjectType::New();

for(it.GoToBegin(); !it.IsAtEnd(); ++it){
if(it.Get() > 0){
InputImageType::IndexType index = it.GetIndex();
//std::cout << "Index: " << index << std::endl;
typename TTransformType::InputPointType    p;
maskImage->TransformIndexToPhysicalPoint(index,p);
//std::cout << "Location: " << p << std::endl;
BlobPointType  blobPoint;
blobPoint.SetPosition(p);
points.push_back(blobPoint);
}
}
 */
} // end namespace itk
#endif
