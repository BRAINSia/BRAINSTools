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

#ifndef _itkTensorLinearInterpolateImageFunction_hxx
#define _itkTensorLinearInterpolateImageFunction_hxx
#include "itkTensorLinearInterpolateImageFunction.h"

#include "itkMath.h"

namespace itk
{
/**
 * Define the number of neighbors
 */
template <typename TInputImage, typename TCoordRep>
const unsigned long TensorLinearInterpolateImageFunction<TInputImage, TCoordRep>::m_Neighbors =
  1 << TInputImage::ImageDimension;

/**
 * Constructor
 */
template <typename TInputImage, typename TCoordRep>
TensorLinearInterpolateImageFunction<TInputImage, TCoordRep>::TensorLinearInterpolateImageFunction() = default;

/**
 * PrintSelf
 */
template <typename TInputImage, typename TCoordRep>
void
TensorLinearInterpolateImageFunction<TInputImage, TCoordRep>::PrintSelf(std::ostream & os, Indent indent) const
{
  this->Superclass::PrintSelf(os, indent);
}

/**
 * Evaluate at image index position
 */
template <typename TInputImage, typename TCoordRep>
typename TensorLinearInterpolateImageFunction<TInputImage, TCoordRep>::OutputType
TensorLinearInterpolateImageFunction<TInputImage, TCoordRep>::EvaluateAtContinuousIndex(
  const ContinuousIndexType & index) const
{
  unsigned int dim; // index over dimension

  /**
   * Compute base index = closet index below point
   * Compute distance from point to base index
   */
  signed long baseIndex[ImageDimension];
  double      distance[ImageDimension];

  for (dim = 0; dim < ImageDimension; dim++)
  {
    baseIndex[dim] = static_cast<long>(floor(index[dim]));
    distance[dim] = index[dim] - static_cast<double>(baseIndex[dim]);
  }

  /**
   * Interpolated value is the weight some of each of the surrounding
   * neighbors. The weight for each neighbour is the fraction overlap
   * of the neighbor pixel with respect to a pixel centered on point.
   */
  OutputType output;
  output.Fill(0.0);

  RealType totalOverlap = 0.0;
  for (unsigned int counter = 0; counter < m_Neighbors; counter++)
  {
    double       overlap = 1.0;   // fraction overlap
    unsigned int upper = counter; // each bit indicates upper/lower neighbour
    IndexType    neighIndex;
    // get neighbor index and overlap fraction
    for (dim = 0; dim < ImageDimension; dim++)
    {
      if (upper & 1)
      {
        neighIndex[dim] = baseIndex[dim] + 1;
        overlap *= distance[dim];
      }
      else
      {
        neighIndex[dim] = baseIndex[dim];
        overlap *= 1.0 - distance[dim];
      }

      upper >>= 1;
    }

    // get neighbor value only if overlap is not zero
    if (overlap)
    {
      const PixelType input = this->GetInputImage()->GetPixel(neighIndex);
      for (unsigned int k = 0; k < 6; k++)
      {
        output[k] += overlap * static_cast<RealType>(input[k]);
      }
      totalOverlap += overlap;
    }

    if (totalOverlap == 1.0)
    {
      // finished
      break;
    }
  }

  return output;
}
} // namespace itk

#endif
