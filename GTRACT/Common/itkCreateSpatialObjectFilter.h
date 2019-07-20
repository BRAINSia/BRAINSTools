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

#ifndef __itkCreateSpatialObjectFilter_h
#define __itkCreateSpatialObjectFilter_h

#include "itkObject.h"
#include "itkImage.h"
#include "itkImageToImageFilter.h"
#include "itkIOCommon.h"
#include "itkPermuteAxesImageFilter.h"
#include "itkVersorRigid3DTransform.h"
#include "itkAffineTransform.h"
#include "itkImageFileReader.h"
#include "itkBlobSpatialObject.h"
#include "itkSpatialObjectPoint.h"
#include "itkSceneSpatialObject.h"

#include <map>
#include <string>

namespace itk
{
/** \class CreateSpatialObjectFilter
 * \brief Creates a Spatial object from a binary image. A transform
 * can be applied to the point coordinates.
 *
 * The output of the filter contains the resulting spatial object.
 */

template < typename TInputImage, typename TTransformType, typename TSpatialObject >
class CreateSpatialObjectFilter : public itk::Object
{
public:
  ITK_DISALLOW_COPY_AND_ASSIGN( CreateSpatialObjectFilter );

  /** Standard class type alias. */
  using Self = CreateSpatialObjectFilter;
  using Superclass = itk::Object;
  using Pointer = SmartPointer< Self >;
  using ConstPointer = SmartPointer< const Self >;

  /** Some convenient type alias. */
  using InputImageType = TInputImage;
  using InputImagePointer = typename InputImageType::Pointer;
  using InputImageConstPointer = typename InputImageType::ConstPointer;
  using InputImageRegionType = typename InputImageType::RegionType;
  using InputImageSizeType = typename InputImageType::SizeType;
  using InputImageSpacingType = typename InputImageType::SpacingType;
  using InputImagePointType = typename InputImageType::PointType;
  using InputImagePixelType = typename InputImageType::PixelType;
  using InputImageIndexType = typename InputImageType::IndexType;
  using InputImageDirectionType = typename InputImageType::DirectionType;

  using TransformType = TTransformType;
  using TransformPointer = typename TransformType::Pointer;

  using SpatialObjectType = TSpatialObject;
  using BlobPointType = typename SpatialObjectType::BlobPointType;
  using SpatialObjectTypePointer = typename SpatialObjectType::Pointer;
  using SpatialObjectPointListType = typename SpatialObjectType::PointListType;

  /** Standard New method. */
  itkNewMacro( Self );

  /** Runtime information support. */
  itkTypeMacro( CreateSpatialObjectFilter, itk::Object );

  /* SetInput and GetOutput Macros */
  itkSetObjectMacro( Input, InputImageType );
  itkGetConstObjectMacro( Output, SpatialObjectType );
  itkSetObjectMacro( Transform, TransformType );
  itkGetMacro( Size, InputImageSizeType );
  itkSetMacro( Size, InputImageSizeType );
  itkGetMacro( Spacing, InputImageSpacingType );
  itkSetMacro( Spacing, InputImageSpacingType );
  itkGetMacro( Origin, InputImagePointType );
  itkSetMacro( Origin, InputImagePointType );

  void
  Update();

protected:
  CreateSpatialObjectFilter();
  ~CreateSpatialObjectFilter() {}

private:
  void
  LoadImage();

  void
  ExtractROI();

  // Input and Output Image
  InputImagePointer        m_Input;
  TransformPointer         m_Transform;
  InputImageSizeType       m_Size;
  InputImageSpacingType    m_Spacing;
  InputImagePointType      m_Origin;
  SpatialObjectTypePointer m_Output;
}; // end of class
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#  include "itkCreateSpatialObjectFilter.hxx"
#endif

#endif
