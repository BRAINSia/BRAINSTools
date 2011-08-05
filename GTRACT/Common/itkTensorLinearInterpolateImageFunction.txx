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

#ifndef _itkTensorLinearInterpolateImageFunction_txx
#define _itkTensorLinearInterpolateImageFunction_txx
#include "itkTensorLinearInterpolateImageFunction.h"

#include "vnl/vnl_math.h"

namespace itk
{
/**
 * Define the number of neighbors
 */
template <class TInputImage, class TCoordRep>
const unsigned long
TensorLinearInterpolateImageFunction<TInputImage, TCoordRep>
::m_Neighbors = 1 << TInputImage::ImageDimension;

/**
 * Constructor
 */
template <class TInputImage, class TCoordRep>
TensorLinearInterpolateImageFunction<TInputImage, TCoordRep>
::TensorLinearInterpolateImageFunction()
{
}

/**
 * PrintSelf
 */
template <class TInputImage, class TCoordRep>
void
TensorLinearInterpolateImageFunction<TInputImage, TCoordRep>
::PrintSelf(std::ostream & os, Indent indent) const
{
  this->Superclass::PrintSelf(os, indent);
}

/**
 * Evaluate at image index position
 */
template <class TInputImage, class TCoordRep>
typename TensorLinearInterpolateImageFunction<TInputImage, TCoordRep>
::OutputType
TensorLinearInterpolateImageFunction<TInputImage, TCoordRep>
::EvaluateAtContinuousIndex(
  const ContinuousIndexType & index) const
{
  unsigned int dim;  // index over dimension

  /**
   * Compute base index = closet index below point
   * Compute distance from point to base index
   */
  signed long baseIndex[ImageDimension];
  double      distance[ImageDimension];

  for( dim = 0; dim < ImageDimension; dim++ )
    {
    baseIndex[dim] = (long)floor( index[dim] );
    distance[dim] = index[dim] - double(baseIndex[dim]);
    }

  /**
   * Interpolated value is the weight some of each of the surrounding
   * neighbors. The weight for each neighbour is the fraction overlap
   * of the neighbor pixel with respect to a pixel centered on point.
   */
  OutputType output;
  output.Fill( 0.0 );

  RealType totalOverlap = 0.0;
  for( unsigned int counter = 0; counter < m_Neighbors; counter++ )
    {
    double       overlap = 1.0;    // fraction overlap
    unsigned int upper = counter;  // each bit indicates upper/lower neighbour
    IndexType    neighIndex;
    // get neighbor index and overlap fraction
    for( dim = 0; dim < ImageDimension; dim++ )
      {
      if( upper & 1 )
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
    if( overlap )
      {
      const PixelType input = this->GetInputImage()->GetPixel( neighIndex );
      for( unsigned int k = 0; k < 6; k++ )
        {
        output[k] += overlap * static_cast<RealType>( input[k] );
        }
      totalOverlap += overlap;
      }

    if( totalOverlap == 1.0 )
      {
      // finished
      break;
      }
    }

  return output;
}
} // namespace itk

#endif
