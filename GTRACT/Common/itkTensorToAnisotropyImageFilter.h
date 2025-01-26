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

#ifndef __itkTensorToAnisotropyImageFilter_h
#define __itkTensorToAnisotropyImageFilter_h

#include "itkObject.h"
#include "itkImage.h"
#include "itkImageToImageFilter.h"
#include "itkIOCommon.h"
#include "gtractCommonWin32.h"

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

enum ENUM_ANISOTROPY_TYPE
{
  MEAN_DIFFUSIVITY = 0,
  FRACTIONAL_ANISOTROPY = 1,
  RELATIVE_ANISOTROPY = 2,
  VOLUME_RATIO = 3,
  AXIAL_DIFFUSIVITY = 4,
  RADIAL_DIFFUSIVITY = 5,
  COHERENCE_INDEX = 6,
  LATTICE_INDEX = 7
};
using AnisotropyType = enum ENUM_ANISOTROPY_TYPE;

class GTRACT_COMMON_EXPORT TensorToAnisotropyImageFilter : public itk::Object
{
public:
  ITK_DISALLOW_COPY_AND_MOVE(TensorToAnisotropyImageFilter);

  /** Standard class type alias. */
  using Self = TensorToAnisotropyImageFilter;
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

  using OutputImageType = itk::Image<float, 3>;
  using OutputImagePointer = OutputImageType::Pointer;
  using OutputImageConstPointer = OutputImageType::ConstPointer;
  using OutputImageRegionType = OutputImageType::RegionType;
  using OutputImagePixelType = OutputImageType::PixelType;
  using OutputImageIndexType = OutputImageType::IndexType;

  /** ImageDimension constants * /
  static constexpr unsigned int InputImageDimension = TInputImage::ImageDimension;
  static constexpr unsigned int OutputImageDimension = TOutputImage::ImageDimension;

  / ** The dimensions of the input image must equal those of the
      output image. * /
  itkConceptMacro(SameDimension,
    (Concept::SameDimension<Self::InputImageDimension,Self::OutputImageDimension>));

  / ** The dimension of the input image must be 4. * /
  itkConceptMacro(DimensionShouldBe4,
    (Concept::SameDimension<Self::InputImageDimension,4>));
*/
  /** Standard New method. */
  itkNewMacro(Self);

  /** Runtime information support. */
  itkOverrideGetNameOfClassMacro(TensorToAnisotropyImageFilter);

  /* SetInput and GetOutput Macros */
  itkSetObjectMacro(Input, InputImageType);
  itkGetConstObjectMacro(Output, OutputImageType);

  itkSetMacro(AnisotropyType, AnisotropyType);

  void
  Update();

protected:
  TensorToAnisotropyImageFilter();
  ~TensorToAnisotropyImageFilter() override = default;

private:
  void
  computVoxelIsotropy();

  void
  computSimpleVoxelAnisotropy();

  void
  computNeighborhoodVoxelAnisotropy();

  // Input and Output Image
  InputImagePointer  m_Input;
  OutputImagePointer m_Output;

  AnisotropyType m_AnisotropyType;
}; // end of class
} // end namespace itk

#endif
