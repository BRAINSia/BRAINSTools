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

#ifndef __itkComputeDiffusionTensorImageFilter_h
#define __itkComputeDiffusionTensorImageFilter_h

#include "itkObject.h"
#include "itkImage.h"
#include "itkImageToImageFilter.h"
#include "itkIOCommon.h"
#include "gtractCommonWin32.h"
#include "algo.h"

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

class GTRACT_COMMON_EXPORT ComputeDiffusionTensorImageFilter : public itk::Object
{
public:
  ITK_DISALLOW_COPY_AND_MOVE(ComputeDiffusionTensorImageFilter);

  /** Standard class type alias. */
  using Self = ComputeDiffusionTensorImageFilter;
  using Superclass = itk::Object;
  using Pointer = SmartPointer<Self>;
  using ConstPointer = SmartPointer<const Self>;

  /** Some convenient type alias. */

  /** Input Image Types */
  using InputImageType = itk::Image<signed short, 4>;
  using InputImagePointer = InputImageType::Pointer;
  using InputImageConstPointer = InputImageType::ConstPointer;
  using InputImageRegionType = InputImageType::RegionType;
  using InputImageSizeType = InputImageType::SizeType;
  using InputImageSpacingType = InputImageType::SpacingType;
  using InputImagePointType = InputImageType::PointType;
  using InputImagePixelType = InputImageType::PixelType;
  using InputImageDirectionType = InputImageType::DirectionType;
  using InputImageIndexType = InputImageType::IndexType;

  /** Output Image Types */
  using OutputPixelType = itk::Vector<float, 6>;
  using OutputImageType = itk::Image<OutputPixelType, 3>;
  using OutputImagePointer = OutputImageType::Pointer;
  using OutputImageConstPointer = OutputImageType::ConstPointer;
  using OutputImageRegionType = OutputImageType::RegionType;
  using OutputImagePixelType = OutputImageType::PixelType;
  using OutputImageIndexType = OutputImageType::IndexType;
  using OutputImageSizeType = OutputImageType::SizeType;
  using OutputImageSpacingType = OutputImageType::SpacingType;
  using OutputImagePointType = OutputImageType::PointType;

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
  itkOverrideGetNameOfClassMacro(ComputeDiffusionTensorImageFilter);

  /* SetInput and GetOutput Macros */
  itkSetObjectMacro(Input, InputImageType);
  itkGetConstObjectMacro(Output, OutputImageType);

  itkSetMacro(UseMedianFilter, bool);
  itkSetMacro(MedianFilterSize, InputImageSizeType);
  itkSetMacro(BackgroundThreshold, int);
  itkSetMacro(NumberOfDirections, int);
  itkSetMacro(NumberOfBSteps, int);
  itkSetMacro(DiffusionDirections, TMatrix);
  itkSetMacro(BValues, TVector);

  itkGetMacro(UseMedianFilter, bool);
  itkGetMacro(MedianFilterSize, InputImageSizeType);
  itkGetMacro(BackgroundThreshold, int);
  itkGetMacro(NumberOfDirections, int);
  itkGetMacro(NumberOfBSteps, int);
  itkGetMacro(DiffusionDirections, TMatrix);
  itkGetMacro(BValues, TVector);

  void
  Update();

protected:
  ComputeDiffusionTensorImageFilter();
  ~ComputeDiffusionTensorImageFilter() override = default;

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

  InputImagePointer m_InternalImage;

  bool               m_UseMedianFilter;
  InputImageSizeType m_MedianFilterSize;
  int                m_BackgroundThreshold;
  int                m_NumberOfDirections;
  int                m_NumberOfBSteps;

  TMatrix m_DiffusionDirections;
  TVector m_BValues;
}; // end of class
} // end namespace itk

#endif
