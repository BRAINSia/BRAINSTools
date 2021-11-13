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
#ifndef itkMixtureStatisticsCostFunction_h_
#define itkMixtureStatisticsCostFunction_h_

#include "itkImageBase.h"
#include "itkImage.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionConstIterator.h"
#include "itkImageRegionConstIteratorWithIndex.h"
#include "itkNumericTraits.h"
#include "itkProgressReporter.h"
#include "itkSingleValuedCostFunction.h"


namespace itk
{
template <typename TFirstImage, typename TSecondImage>
class MixtureStatisticCostFunction : public SingleValuedCostFunction
{
public:
  ITK_DISALLOW_COPY_AND_MOVE(MixtureStatisticCostFunction);

  /** Standard type alias. */
  using Self = MixtureStatisticCostFunction;
  using Superclass = SingleValuedCostFunction;
  using Pointer = SmartPointer<Self>;
  using ConstPointer = SmartPointer<const Self>;

  /** Run-time type information (and related methods).   */
  itkTypeMacro(MixtureStatisticCostFunction, MultipleValuedCostFunction);

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /**  Type of the First Image. */
  using FirstImageType = TFirstImage;
  using FirstImageConstPointer = typename FirstImageType::ConstPointer;

  /**  Type of the Second Image. */
  using SecondImageType = TSecondImage;
  using SecondImageConstPointer = typename SecondImageType::ConstPointer;

  /** Array Typedefs. */
  using ParametersType = Superclass::ParametersType;
  using MeasureType = Superclass::MeasureType;
  using DerivativeType = Superclass::DerivativeType;

  /**  Type for the mask of the first image. Only pixels that are "inside"
   * this mask will be considered for the computation of the metric */
  using ImageMaskType = itk::Image<signed short, 3>;
  using ImageMaskPointer = typename ImageMaskType::Pointer;

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

  /** The dimensions of parameter space. */
  static constexpr unsigned int number_of_unknowns = 1;

  /** Not necessary for this optimizer. */
  void
  GetDerivative(const ParametersType & itkNotUsed(parameters), DerivativeType & itkNotUsed(derivative)) const override
  {}

  /** Return the values evaluated for the given parameters. */
  MeasureType
  GetValue(const ParametersType & parameters) const override;

  /** Get the SpaceDimension. */
  unsigned int
  GetNumberOfParameters() const override;

protected:
  MixtureStatisticCostFunction();
  ~MixtureStatisticCostFunction() override;
  void
  PrintSelf(std::ostream & os, Indent indent) const override;

  FirstImageConstPointer  m_FirstImage;
  SecondImageConstPointer m_SecondImage;

  mutable ImageMaskPointer m_ImageMask;

private:
  /** Measurement value. */
  mutable MeasureType m_Measure;
};
} // end namespace itk

#include "itkMixtureStatisticCostFunction.hxx"
#endif // itkMixtureStatisticsCostFunction_h_
