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
typedef enum ENUM_ANISOTROPY_TYPE AnisotropyType;

class GTRACT_COMMON_EXPORT TensorToAnisotropyImageFilter : public itk::Object
{
public:
  /** Standard class typedefs. */
  typedef TensorToAnisotropyImageFilter Self;
  typedef itk::Object                   Superclass;
  typedef SmartPointer<Self>            Pointer;
  typedef SmartPointer<const Self>      ConstPointer;

  /** Some convenient typedefs. */
  typedef itk::Vector<float, 6>         InputPixelType;
  typedef itk::Image<InputPixelType, 3> InputImageType;
  typedef InputImageType::Pointer       InputImagePointer;
  typedef InputImageType::ConstPointer  InputImageConstPointer;
  typedef InputImageType::RegionType    InputImageRegionType;
  typedef InputImageType::SizeType      InputImageSizeType;
  typedef InputImageType::SpacingType   InputImageSpacingType;
  typedef InputImageType::PointType     InputImagePointType;
  typedef InputImageType::PixelType     InputImagePixelType;
  typedef InputImageType::DirectionType InputImageDirectionType;

  typedef itk::Image<float, 3>          OutputImageType;
  typedef OutputImageType::Pointer      OutputImagePointer;
  typedef OutputImageType::ConstPointer OutputImageConstPointer;
  typedef OutputImageType::RegionType   OutputImageRegionType;
  typedef OutputImageType::PixelType    OutputImagePixelType;
  typedef OutputImageType::IndexType    OutputImageIndexType;

  /** ImageDimension constants * /
  itkStaticConstMacro(InputImageDimension, unsigned int,
                      TInputImage::ImageDimension);
  itkStaticConstMacro(OutputImageDimension, unsigned int,
                      TOutputImage::ImageDimension);

  / ** The dimensions of the input image must equal those of the
      output image. * /
  itkConceptMacro(SameDimension,
    (Concept::SameDimension<itkGetStaticConstMacro(InputImageDimension),itkGetStaticConstMacro(OutputImageDimension)>));

  / ** The dimension of the input image must be 4. * /
  itkConceptMacro(DimensionShouldBe4,
    (Concept::SameDimension<itkGetStaticConstMacro(InputImageDimension),4>));
*/
  /** Standard New method. */
  itkNewMacro(Self);

  /** Runtime information support. */
  itkTypeMacro(TensorToAnisotropyImageFilter, itk::Object);

  /* SetInput and GetOutput Macros */
  itkSetObjectMacro(Input,  InputImageType);
  itkGetConstObjectMacro(Output, OutputImageType);

  itkSetMacro(AnisotropyType, AnisotropyType);

  void Update();

protected:
  TensorToAnisotropyImageFilter();
  ~TensorToAnisotropyImageFilter()
  {
  }

private:
  ITK_DISALLOW_COPY_AND_ASSIGN(TensorToAnisotropyImageFilter);

  void computVoxelIsotropy();

  void computSimpleVoxelAnisotropy();

  void computNeighborhoodVoxelAnisotropy();

  // Input and Output Image
  InputImagePointer  m_Input;
  OutputImagePointer m_Output;

  AnisotropyType m_AnisotropyType;
};  // end of class
} // end namespace itk

#endif
