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
  /** Standard class typedefs. */
  typedef ComputeDiffusionTensorImageFilter Self;
  typedef itk::Object                       Superclass;
  typedef SmartPointer<Self>                Pointer;
  typedef SmartPointer<const Self>          ConstPointer;

  /** Some convenient typedefs. */

  /** Input Image Types */
  typedef itk::Image<signed short, 4>   InputImageType;
  typedef InputImageType::Pointer       InputImagePointer;
  typedef InputImageType::ConstPointer  InputImageConstPointer;
  typedef InputImageType::RegionType    InputImageRegionType;
  typedef InputImageType::SizeType      InputImageSizeType;
  typedef InputImageType::SpacingType   InputImageSpacingType;
  typedef InputImageType::PointType     InputImagePointType;
  typedef InputImageType::PixelType     InputImagePixelType;
  typedef InputImageType::DirectionType InputImageDirectionType;
  typedef InputImageType::IndexType     InputImageIndexType;

  /** Output Image Types */
  typedef itk::Vector<float, 6>          OutputPixelType;
  typedef itk::Image<OutputPixelType, 3> OutputImageType;
  typedef OutputImageType::Pointer       OutputImagePointer;
  typedef OutputImageType::ConstPointer  OutputImageConstPointer;
  typedef OutputImageType::RegionType    OutputImageRegionType;
  typedef OutputImageType::PixelType     OutputImagePixelType;
  typedef OutputImageType::IndexType     OutputImageIndexType;
  typedef OutputImageType::SizeType      OutputImageSizeType;
  typedef OutputImageType::SpacingType   OutputImageSpacingType;
  typedef OutputImageType::PointType     OutputImagePointType;

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
  itkTypeMacro(ComputeDiffusionTensorImageFilter, itk::Object);

  /* SetInput and GetOutput Macros */
  itkSetObjectMacro(Input,  InputImageType);
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

  void Update();

protected:
  ComputeDiffusionTensorImageFilter();
  ~ComputeDiffusionTensorImageFilter()
  {
  }

private:
  ITK_DISALLOW_COPY_AND_ASSIGN(ComputeDiffusionTensorImageFilter);

  void computVoxelIsotropy();

  void computSimpleVoxelAnisotropy();

  void computNeighborhoodVoxelAnisotropy();

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
};  // end of class
} // end namespace itk

#endif
