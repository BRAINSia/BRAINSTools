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

#ifndef __itkInvertBSplineFilter_h
#define __itkInvertBSplineFilter_h

#include "itkObject.h"
#include "itkImage.h"
#include "itkImageToImageFilter.h"
#include "itkIOCommon.h"
#include <itkExtractImageFilter.h>
#include <itkMetaDataDictionary.h>
#include <itkMetaDataObject.h>
#include <itkImageRegistrationMethod.h>
#include <itkMattesMutualInformationImageToImageMetric.h>
#include <itkLinearInterpolateImageFunction.h>
#include <itkBSplineTransform.h>
#include <itkThinPlateR2LogRSplineKernelTransform.h>
#include <itkLBFGSBOptimizer.h>
#include <itkCenteredTransformInitializer.h>
#include <itkTimeProbesCollectorBase.h>
#include <itkTransformFactory.h>
#include "gtractCommonWin32.h"

#include <map>
#include <string>

namespace itk
{
/** \class InvertBSplineFilter
 *
 */

class GTRACT_COMMON_EXPORT InvertBSplineFilter : public itk::Object
{
public:
  /** Standard class typedefs. */
  typedef InvertBSplineFilter      Self;
  typedef itk::Object              Superclass;
  typedef SmartPointer<Self>       Pointer;
  typedef SmartPointer<const Self> ConstPointer;

  static const unsigned int transformDimension = 3;
  typedef double CoordinateRepresentationType;

  /** Fixed Image typedefs. */
  typedef itk::Image<signed short, transformDimension> ImageType;
  typedef ImageType::Pointer                           ImageTypePointer;
  typedef ImageType::ConstPointer                      ImageConstPointer;
  typedef ImageType::RegionType                        ImageRegionType;
  typedef ImageType::SizeType                          ImageSizeType;
  typedef ImageType::SpacingType                       ImageSpacingType;
  typedef ImageType::PointType                         ImagePointType;
  typedef ImageType::PixelType                         ImagePixelType;
  typedef ImageType::DirectionType                     ImageDirectionType;
  typedef ImageType::IndexType                         ImageIndexType;

  /** B-Spline Transform typedefs */
  static const unsigned int SplineOrder = 3;
  typedef itk::BSplineTransform<
      CoordinateRepresentationType,
      transformDimension,
      SplineOrder>     BsplineTransformType;
  typedef BsplineTransformType::Pointer BsplineTransformTypePointer;

  /** Output Transform typedefs. */
  typedef itk::ThinPlateR2LogRSplineKernelTransform<CoordinateRepresentationType,
                                                    transformDimension>  TransformType;
  typedef TransformType::Pointer                                       TransformTypePointer;
  typedef itk::Point<CoordinateRepresentationType, transformDimension> PointType;
  typedef TransformType::PointSetType                                  PointSetType;
  typedef PointSetType::PointIdentifier                                PointIdType;

  /** Standard New method. */
  itkNewMacro(Self);

  /** Runtime information support. */
  itkTypeMacro(InvertBSplineFilter, itk::Object);

  /* SetInput and GetOutput Macros */
  itkSetObjectMacro(Input, BsplineTransformType);
  itkSetObjectMacro(ExampleImage, ImageType);
  itkGetConstObjectMacro(Output, TransformType);

  itkSetMacro(XgridSize, int);
  itkSetMacro(YgridSize, int);
  itkSetMacro(ZgridSize, int);

  itkGetMacro(XgridSize, int);
  itkGetMacro(YgridSize, int);
  itkGetMacro(ZgridSize, int);

  void Update();

protected:
  InvertBSplineFilter();
  ~InvertBSplineFilter()
  {
  }

private:
  ITK_DISALLOW_COPY_AND_ASSIGN(InvertBSplineFilter);

  /*** Input and Output Objects ***/
  BsplineTransformTypePointer m_Input;
  TransformTypePointer        m_Output;
  ImageTypePointer            m_ExampleImage;

  int m_XgridSize;
  int m_YgridSize;
  int m_ZgridSize;
};      // end of class
} // end namespace itk

#endif
