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

#ifndef __itkFastMarchingCostFunction_h
#define __itkFastMarchingCostFunction_h

#include "itkSingleValuedCostFunction.h"

#include <itkMath.h>
#include <iostream>
#include <fstream>
#include <itkImage.h>
#include <itkObject.h>
#include <itkImageToImageFilter.h>
#include <itkIOCommon.h>
#include "itkProcessObject.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkImageRegionConstIteratorWithIndex.h"
#include "GtractTypes.h"

#include <map>
#include <string>

#include "itkIndex.h"
#include "itkMath.h"

namespace itk
{
class FastMarchingCostFunction : public SingleValuedCostFunction
{
public:
  ITK_DISALLOW_COPY_AND_ASSIGN(FastMarchingCostFunction);

  /** Standard class type alias */
  using Self = FastMarchingCostFunction;
  using Superclass = SingleValuedCostFunction;
  using Pointer = SmartPointer<Self>;
  using ConstPointer = SmartPointer<const Self>;

  /** Run-time type information (and related methods). */
  itkTypeMacro(FastMarchingCostFunction, SingleValuedCostFunction);

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Some convenient type alias. */
  using CostImageType = itk::Image<float, 3>;
  using CostImagePointer = CostImageType::Pointer;
  using CostImageConstPointer = CostImageType::ConstPointer;
  using CostImageRegionType = CostImageType::RegionType;
  using CostImageSizeType = CostImageType::SizeType;
  using CostImageSpacingType = CostImageType::SpacingType;
  using CostImagePointType = CostImageType::PointType;
  using CostImagePixelType = CostImageType::PixelType;
  using CostImageIndexType = CostImageType::IndexType;
  using CostImageDirectionType = CostImageType::DirectionType;

  using CostIPType = itk::LinearInterpolateImageFunction<CostImageType, float>;      //
                                                                                     //
                                                                                     // ScalarIPType;
  using CostIPTypePointer = CostIPType::Pointer;

  /** ImageDimension constants */
  static constexpr unsigned int CostImageDimension = 3;

  /*  A position in the optimization space. */
  using ParametersType = Superclass::ParametersType;

  /*  std::cost function derivative (gradient). */
  using DerivativeType = Superclass::DerivativeType;

  /*  std::cost function value. */
  using MeasureType = Superclass::MeasureType;

  /** Set/Get input Cost Image  */
  itkSetObjectMacro( CostImage, CostImageType );
  itkGetConstObjectMacro( CostImage, CostImageType );

  // Returns dimension of image
  unsigned int GetNumberOfParameters() const override;

  /** This method returns the value of the std::cost function for
    * the specified parameters, or position. */
  MeasureType GetValue( const ParametersType & parameters ) const override;

  /** This method returns the derivative of the std::cost function corresponding
    * to the specified parameters.   */
  void GetDerivative( const ParametersType & parameters, DerivativeType & derivative ) const override;

protected:
  FastMarchingCostFunction();
  // virtual ~FastMarchingCostFunction(){};
  ~FastMarchingCostFunction() override
  {
  }

  void PrintSelf(std::ostream & os, Indent indent) const override;

  // void SetMetaDataHeader();

  CostImagePointer    m_CostImage;
  CostImageRegionType m_BufferedRegion;
  CostIPTypePointer   m_CostIP;
private:

};                                        // end of class
} // end of namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkFastMarchingCostFunction.cxx"
#endif

#endif
