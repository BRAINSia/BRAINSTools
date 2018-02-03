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

#ifndef _itkTensorInterpolateImageFunction_h
#define _itkTensorInterpolateImageFunction_h
#include "itkFixedArray.h"
#include "itkImageFunction.h"
#include "itkSymmetricSecondRankTensor.h"

namespace itk
{
template <
  class TInputImage,
  class TCoordRep = float // ,typename TComponent = float
  >
class TensorInterpolateImageFunction :
  public         ImageFunction<
    TInputImage,
    SymmetricSecondRankTensor<double, 3>,
    TCoordRep>
{
public:
  static constexpr unsigned int ImageDimension = TInputImage::ImageDimension;

  /** Standard class typedefs. */
  typedef TensorInterpolateImageFunction Self;
  typedef ImageFunction<TInputImage,
                        SymmetricSecondRankTensor<double, 3>, TCoordRep> Superclass;
  typedef SmartPointer<Self>       Pointer;
  typedef SmartPointer<const Self> ConstPointer;

  /** Run-time type information (and related methods). */
  itkTypeMacro(TensorInterpolateImageFunction, ImageFunction);

  /** InputImageType typedef support. */
  typedef typename Superclass::InputImageType                  InputImageType;
  typedef typename InputImageType::PixelType                   PixelType;
  typedef itk::SymmetricSecondRankTensor<double, 3>::ValueType ValueType;

  // typedef typename PixelType::ValueType       ValueType;
  typedef typename NumericTraits<ValueType>::RealType RealType;

  /** Point typedef support. */
  typedef typename Superclass::PointType PointType;

  /** Index typedef support. */
  typedef typename Superclass::IndexType IndexType;

  /** ContinuousIndex typedef support. */
  typedef typename Superclass::ContinuousIndexType ContinuousIndexType;

  /** Output type is SymmetricSecondRankTensor<TComponent,Dimension>. */
  typedef typename Superclass::OutputType OutputType;

  /** CoordRep typedef support. */
  typedef TCoordRep CoordRepType;
  OutputType Evaluate( const PointType & point ) const override
  {
    ContinuousIndexType index;

    this->GetInputImage()->TransformPhysicalPointToContinuousIndex( point, index );
    return this->EvaluateAtContinuousIndex( index );
  }

  /** Interpolate the image at a continuous index position
   *
   * Returns the interpolated image intensity at a
   * specified index position. No bounds checking is done.
   * The point is assume to lie within the image buffer.
   *
   * Subclasses must override this method.
   *
   * ImageFunction::IsInsideBuffer() can be used to check bounds before
   * calling the method. */
  OutputType EvaluateAtContinuousIndex(constexpr ContinuousIndexType & index ) constexpr override   = 0;

  /** Interpolate the image at an index position.
   * Simply returns the image value at the
   * specified index position. No bounds checking is done.
   * The point is assume to lie within the image buffer.
   *
   * ImageFunction::IsInsideBuffer() can be used to check bounds before
   * calling the method. */
  OutputType EvaluateAtIndex( const IndexType & index ) const override
  {
    OutputType output;
    PixelType  input = this->GetInputImage()->GetPixel( index );

    for( unsigned int k = 0; k < 6; k++ )
      {
      output[k] = static_cast<double>( input[k] );
      }
    return output;
  }

protected:
  TensorInterpolateImageFunction()
  {
  }

  ~TensorInterpolateImageFunction() override
  {
  }

  void PrintSelf(std::ostream & os, Indent indent) const override
  {
    Superclass::PrintSelf( os, indent );
  }

private:
  ITK_DISALLOW_COPY_AND_ASSIGN(TensorInterpolateImageFunction);
};
} // namespace itk

#endif
