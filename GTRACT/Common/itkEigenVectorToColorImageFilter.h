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

#ifndef __itkEigenVectorToColorImageFilter_h
#define __itkEigenVectorToColorImageFilter_h

#include "itkObject.h"
#include "itkImage.h"
#include "itkImageToImageFilter.h"
#include "itkIOCommon.h"
#include "itkRGBAPixel.h"

#include <map>
#include <string>

namespace itk
{
/** \class TensorToAnisotropyImageFilter
 * \brief Calculates the Specified Anisotropy Index.
 *
 * The following Anisotropy Image are supported:
 *    Fractional Anistropy
 *    Relatibve Anisotropy
 *    Volume Ratio
 *    Radial Diffusivity
 *    Axial Diffusivity
 *    Coheernce Index
 *    Lattice Index
 *    Mean Diffusivity
 *
 */

enum ENUM_TENSOR_SHAPE_TYPE
  {
  TERTIARY_EIGENVECTOR = 0,
  SECONDARY_EIGENVECTOR = 1,
  PRIMARY_EIGENVECTOR = 2,
  TENSOR_SHAPE = 3
  };
using TensorShapeType = enum ENUM_TENSOR_SHAPE_TYPE;

class EigenVectorToColorImageFilter : public itk::Object
{
public:
  ITK_DISALLOW_COPY_AND_ASSIGN(EigenVectorToColorImageFilter);

  /** Standard class type alias. */
  using Self = EigenVectorToColorImageFilter;
  using Superclass = itk::Object;
  using Pointer = SmartPointer<Self>;
  using ConstPointer = SmartPointer<const Self>;

  /** Some convenient type alias. */
  using InputPixelType = itk::Vector<float, 6>;
  using InputImageType = itk::Image<InputPixelType, 3>;
  using InputImagePointer = InputImageType::Pointer;
  using InputImageConstPointer = InputImageType::ConstPointer;
  using InputImageRegionType = InputImageType::RegionType;
  using InputImageSizeType = InputImageType::SizeType;
  using InputImageSpacingType = InputImageType::SpacingType;
  using InputImagePointType = InputImageType::PointType;
  using InputImagePixelType = InputImageType::PixelType;
  using InputImageDirectionType = InputImageType::DirectionType;

  using OutputPixelType = itk::RGBAPixel<unsigned char>;
  using OutputImageType = itk::Image<OutputPixelType, 3>;
  using OutputImagePointer = OutputImageType::Pointer;
  using OutputImageRegionType = OutputImageType::RegionType;


/** Standard New method. */
  itkNewMacro(Self);

/** Runtime information support. */
  itkTypeMacro(EigenVectorToColorImageFilter, itk::Object);

  itkSetObjectMacro(Input,  InputImageType);
  itkGetConstObjectMacro(Output, OutputImageType);
  itkSetMacro(TensorShapeType, TensorShapeType);

  void Update();

protected:
  EigenVectorToColorImageFilter();
  ~EigenVectorToColorImageFilter() override
  {
  }

private:
// Input and Output Image
  InputImagePointer  m_Input;
  OutputImagePointer m_Output;

  TensorShapeType m_TensorShapeType;
};  // end of class
} // end namespace itk

#endif
