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
template <class TInputImage, class TTransformType, class TSpatialObject>
CreateSpatialObjectFilter<TInputImage, TTransformType, TSpatialObject>::CreateSpatialObjectFilter()
{
  m_Input = NULL;
  m_Output = NULL;
  m_Transform = NULL;
}

template <class TInputImage, class TTransformType, class TSpatialObject>
void CreateSpatialObjectFilter<TInputImage, TTransformType, TSpatialObject>::Update()
{
  SpatialObjectPointListType points;

  m_Output  = SpatialObjectType::New();

  typedef itk::ImageRegionConstIteratorWithIndex<InputImageType> IteratorType;
  IteratorType it( m_Input, m_Input->GetLargestPossibleRegion() );
  // int count;
  for( it.GoToBegin(); !it.IsAtEnd(); ++it )
    {
    if( it.Get() > 0 )
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
typedef itk::ResampleImageFilter<
InputImageType,
InputImageType >    ResampleFilterType;

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



typedef itk::ImageRegionConstIteratorWithIndex< InputImageType >  IteratorType;
IteratorType it(maskImage, maskImage->GetLargestPossibleRegion());

typedef itk::SpatialObjectPoint<3>    BlobPointType;
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
