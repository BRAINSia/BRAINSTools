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
 *
 *  Program:   Insight Segmentation & Registration Toolkit
 *  Module:    $RCSfile: itkMixtureStatisticsCostFunction.h,v $
 *  Language:  C++
 *  Date:      $Date: 2007-03-29 19:37:00 $
 *  Version:   $Revision: 1.14 $
 *
 *  Copyright (c) Insight Software Consortium. All rights reserved.
 *  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.
 *
 *  This software is distributed WITHOUT ANY WARRANTY; without even
 *  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 *  PURPOSE.  See the above copyright notices for more information.
 *
 *  =========================================================================*/
#ifndef __itkMixtureStatisticsCostFunction_h
#define __itkMixtureStatisticsCostFunction_h

#include "itkImageBase.h"
#include "itkImage.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionConstIterator.h"
#include "itkImageRegionConstIteratorWithIndex.h"
#include "itkNumericTraits.h"
#include "itkProgressReporter.h"
#include "itkMultipleValuedCostFunction.h"
#include "itkExceptionObject.h"

namespace itk
{
template <typename TFirstImage, typename TSecondImage>
class MixtureStatisticCostFunction :
  public MultipleValuedCostFunction
{
public:
  ITK_DISALLOW_COPY_AND_ASSIGN(MixtureStatisticCostFunction);

  /** Standard type alias. */
  using Self = MixtureStatisticCostFunction;
  using Superclass = MultipleValuedCostFunction;
  using Pointer = SmartPointer<Self>;
  using ConstPointer = SmartPointer<const Self>;

  /** Run-time type information (and related methods).   */
  itkTypeMacro(MixtureStatisticCostFunction, MultipleValuedCostFunction);

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /**  Type of the First Image. */
  using FirstImageType = TFirstImage;
  using FirstImagePixelType = typename TFirstImage::PixelType;
  using FirstImageConstPointer = typename FirstImageType::ConstPointer;

  /**  Type of the Second Image. */
  using SecondImageType = TSecondImage;
  using SecondImagePixelType = typename TFirstImage::PixelType;
  using SecondImageConstPointer = typename SecondImageType::ConstPointer;
  using SecondImageRegionType = typename SecondImageType::RegionType;

  /** Constants for the image dimensions */
  static constexpr unsigned int FirstImageDimension = TFirstImage::ImageDimension;
  static constexpr unsigned int SecondImageDimension = TSecondImage::ImageDimension;

  /** Array Typedefs. */
  using ParametersType = Superclass::ParametersType;
  using MeasureType = Superclass::MeasureType;
  using DerivativeType = Superclass::DerivativeType;

  /**  Type for the mask of the first image. Only pixels that are "inside"
    * this mask will be considered for the computation of the metric */
  using ImageMaskType = itk::Image<signed short, 3>;
  using ImageMaskPointer = typename  ImageMaskType::Pointer;

  /** Connect the First Image.  */
  itkSetConstObjectMacro(FirstImage, FirstImageType);

  /** Get the First Image. */
  itkGetConstObjectMacro(FirstImage, FirstImageType);

  /** Connect the Second Image.  */
  itkSetConstObjectMacro(SecondImage, SecondImageType);

  /** Get the Second Image. */
  itkGetConstObjectMacro(SecondImage, SecondImageType);

  /** Set/Get the first image mask. */
  itkSetObjectMacro(ImageMask, ImageMaskType);
  itkGetConstObjectMacro(ImageMask, ImageMaskType);

  itkGetMacro(DesiredMean, double);
  itkSetMacro(DesiredMean, double);
  itkGetMacro(DesiredVariance, double);
  itkSetMacro(DesiredVariance, double);

  itkGetMacro(NumberOfMaskVoxels, double);
  itkGetMacro(SumOfFirstMaskVoxels, double);
  itkGetMacro(SumOfSecondMaskVoxels, double);
  itkGetMacro(SumSquaresOfFirstMaskVoxels, double);
  itkGetMacro(SumSquaresOfSecondMaskVoxels, double);
  itkGetMacro(SumOfFirstTimesSecondMaskVoxels, double);

  /** The dimensions of parameter space. */
  enum { SpaceDimension = 2 };

  /** Not necessary for this optimizer. */
  void GetDerivative( const ParametersType & itkNotUsed(parameters),
                      DerivativeType & itkNotUsed(derivative) ) const override
  {
  }

  /** Return the values evaluated for the given parameters. */
  MeasureType GetValue(const ParametersType & parameters) const override;

  /** Return a pointer of values evaluated for the given parameters. */
  MeasureType * GetValue(ParametersType & parameters);

  /** Get the SpaceDimension. */
  unsigned int GetNumberOfParameters() const override;

  /** Get the number Range Dimension. */
  unsigned int GetNumberOfValues() const override;

  /** Initialize */
  void Initialize(short label);

protected:
  MixtureStatisticCostFunction();
  ~MixtureStatisticCostFunction() override;

  void PrintSelf(std::ostream & os, Indent indent) const override;

  FirstImageConstPointer  m_FirstImage;
  SecondImageConstPointer m_SecondImage;

  mutable ImageMaskPointer m_ImageMask;
private:
  double m_DesiredMean;
  double m_DesiredVariance;

  double m_NumberOfMaskVoxels;
  double m_SumOfFirstMaskVoxels;
  double m_SumOfSecondMaskVoxels;
  double m_SumSquaresOfFirstMaskVoxels;
  double m_SumSquaresOfSecondMaskVoxels;
  double m_SumOfFirstTimesSecondMaskVoxels;

  /** Different arrays. */
  mutable MeasureType    m_Measure;
  mutable MeasureType *  m_MeasurePointer;
  mutable ParametersType m_Parameters;
};
}   // end namespace itk
#ifndef ITK_MANUAL_INSTANTIATION
#include "itkMixtureStatisticCostFunction.hxx"
#endif

#endif
