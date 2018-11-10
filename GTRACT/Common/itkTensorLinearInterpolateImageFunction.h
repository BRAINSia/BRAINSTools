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

#ifndef _itkTensorLinearInterpolateImageFunction_h
#define _itkTensorLinearInterpolateImageFunction_h

#include "itkTensorInterpolateImageFunction.h"

namespace itk
{
template <typename TInputImage, typename TCoordRep = float>
class TensorLinearInterpolateImageFunction :
  public         TensorInterpolateImageFunction<TInputImage, TCoordRep>
{
public:
  /** Standard class type alias. */
  using Self = TensorLinearInterpolateImageFunction;
  using Superclass = TensorInterpolateImageFunction<TInputImage, TCoordRep>;
  using Pointer = SmartPointer<Self>;
  using ConstPointer = SmartPointer<const Self>;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(TensorLinearInterpolateImageFunction,
               TensorInterpolateImageFunction);

  /** InputImageType type alias support. */
  using InputImageType = typename Superclass::InputImageType;
  using PixelType = typename Superclass::PixelType;
  using ValueType = typename Superclass::ValueType;
  using RealType = typename Superclass::RealType;

  /** Grab the vector dimension from the superclass. */
  // static constexpr unsigned int Dimension = //    Superclass::Dimension;

  /** Dimension underlying input image. */
  static constexpr unsigned int ImageDimension = Superclass::ImageDimension;

  /** Index type alias support. */
  using IndexType = typename Superclass::IndexType;

  /** ContinuousIndex type alias support. */
  using ContinuousIndexType = typename Superclass::ContinuousIndexType;

  /** Output type is Vector<double,Dimension> */
  using OutputType = typename Superclass::OutputType;

  /** Evaluate the function at a ContinuousIndex position
   *
   * Returns the linearly interpolated image intensity at a
   * specified point position. No bounds checking is done.
   * The point is assume to lie within the image buffer.
   *
   * ImageFunction::IsInsideBuffer() can be used to check bounds before
   * calling the method. */
  OutputType EvaluateAtContinuousIndex(const ContinuousIndexType & index ) const override;

protected:
  TensorLinearInterpolateImageFunction();
  ~TensorLinearInterpolateImageFunction() override
  {
  }

  void PrintSelf(std::ostream & os, Indent indent) const override;

private:
  TensorLinearInterpolateImageFunction(const Self &); // purposely not
                                                      // implemented
  void operator=(const Self &);                       // purposely not

  // implemented

  /** Number of neighbors used in the interpolation */
  static const unsigned long m_Neighbors;
};
} // namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkTensorLinearInterpolateImageFunction.hxx"
#endif

#endif
