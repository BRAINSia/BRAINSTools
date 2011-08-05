/*=========================================================================

 Program:   GTRACT (Guided Tensor Restore Anatomical Connectivity Tractography)
 Module:    $RCSfile: $
 Language:  C++
 Date:      $Date: 2006/03/29 14:53:40 $
 Version:   $Revision: 1.9 $

   Copyright (c) University of Iowa Department of Radiology. All rights reserved.
   See GTRACT-Copyright.txt or http://mri.radiology.uiowa.edu/copyright/GTRACT-Copyright.txt
   for details.

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

      This software is distributed WITHOUT ANY WARRANTY; without even
      the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
      PURPOSE.  See the above copyright notices for more information.

=========================================================================*/

#ifndef __itkFastMarchingCostFunction_cxx
#define __itkFastMarchingCostFunction_cxx

#include "itkFastMarchingCostFunction.h"

#include <iostream>
#include <vector>

#include <itkImageRegionIterator.h>
#include <itkImageRegionConstIterator.h>
#include <itkIOCommon.h>

#include <itkImageRegionConstIteratorWithIndex.h>
#include "itkConstNeighborhoodIterator.h"
#include <itkNeighborhoodAlgorithm.h>
#include "itkLinearInterpolateImageFunction.h"
#include <itkVector.h>
#include <itkListSample.h>
#include <vnl/vnl_vector.h>
#include <itkFixedArray.h>
// #include "itkMetaDataObject.h"

#include "itkNumericTraits.h"
#include "vnl/vnl_math.h"
#include <algorithm>

namespace itk
{
// template <class TCostImage >
FastMarchingCostFunction // < TCostImage >
::FastMarchingCostFunction()
{
  m_CostIP = CostIPType::New();
}

unsigned int
FastMarchingCostFunction
::GetNumberOfParameters() const
{
  unsigned int size = 0;

  size = CostImageDimension;
  return size;
  // return SpaceDimension;
}

// template <class TCostImage>
void
FastMarchingCostFunction // < TCostImage>
::PrintSelf(std::ostream & os, Indent indent) const
{
  Superclass::PrintSelf( os, indent );

  os << indent << "Input Cost Image: " << m_CostImage.GetPointer()  << std::endl;
}

// template <class TCostImage>
// typename FastMarchingCostFunction< TCostImage>::MeasureType
FastMarchingCostFunction::MeasureType
FastMarchingCostFunction
::GetValue(const ParametersType & parameters) const
{
  /* Cost Function = |gradient(T)| = |T(r)-T(n)| = |r-n|/F(r):
    where "T" is the vcl_cost, "r" is initial point,"n" is neighbor.
    F(r) is speed of initial point, and |r-n| is distance between the points*/

  float initialCost;

  CostIPType::ContinuousIndexType inputIndex, neighIndex;

  m_CostIP->SetInputImage(m_CostImage);
  CostImageRegionType costRegion = m_CostImage->GetLargestPossibleRegion();
  double              value = 0.0;
  for( unsigned int i = 0; i < CostImageDimension; i++ )
    {
    inputIndex[i] = parameters[i];
    }

  initialCost = m_CostIP->EvaluateAtContinuousIndex(inputIndex);

  value = initialCost;

  return value;
} // end of GetValue

// template <class TCostImage>
void
FastMarchingCostFunction // < TCostImage>
::GetDerivative( const ParametersType & parameters,
                 DerivativeType  & derivative ) const
{
  float initialCost;

  CostIPType::ContinuousIndexType inputIndex, neighIndex;

  m_CostIP->SetInputImage(m_CostImage);
  CostImageRegionType costRegion = m_CostImage->GetLargestPossibleRegion();
  for( unsigned int i = 0; i < CostImageDimension; i++ )
    {
    inputIndex[i] = parameters[i];
    }

  if( !costRegion.IsInside( inputIndex ) )  // not in Region so skip, normal
                                            // will be 0.0
    {
    std::cout << "Warning initial point is outside of region " << std::endl;
    }

  typedef vnl_vector_fixed<float, CostImageDimension> FVector;
  FVector offset;

  FVector sum; sum.fill(0);
  FVector newSum; newSum.fill(0);
  FVector neighOffset;              // vector representation of neighborhood
                                    // offsets
  FVector normal;   normal.fill(0); // normal used for derivative
  double  value = 0.0;              // Cost Function solution

  CostImageSpacingType spacing = m_CostImage->GetSpacing();

  initialCost = m_CostIP->EvaluateAtContinuousIndex(inputIndex);

  value = initialCost;

  // Get complete neighborhood of selected point to calculate normal (direction)

  float neighCost = 0.0;
  float lowNeighCost = initialCost;
  for( int i = -1; i < 2; i++ )
    {
    for( int j = -1; j < 2; j++ )
      {
      for( int k = -1; k < 2; k++ )  // Fix later to account for dimension
        {
        neighIndex[0] = inputIndex[0] + static_cast<float>( i ) / 1; //
                                                                     // +i,+i/1.5,
                                                                     // /4
        neighIndex[1] = inputIndex[1] + static_cast<float>( j ) / 1;
        neighIndex[2] = inputIndex[2] + static_cast<float>( k ) / 1;

        if( neighIndex == inputIndex )  // same index so skip this neighbor
          {
          continue;
          }

        if( !costRegion.IsInside( neighIndex ) )  // not in Region so skip
          {
          continue;
          }
        offset[0] = static_cast<float>( i ) / 1;
        offset[1] = static_cast<float>( j ) / 1;
        offset[2] = static_cast<float>( k ) / 1;
        /* Compute distance from initial point to its neighbor */
        for( unsigned int n = 0; n < CostImageDimension; n++ )
          {
          neighOffset[n] = offset[n] * vnl_math_abs(spacing[n]);
          }

        neighCost = m_CostIP->EvaluateAtContinuousIndex(neighIndex);

        if( neighCost < lowNeighCost )
          {
          lowNeighCost = neighCost;
          newSum = neighOffset;
          }

        if( neighCost < initialCost )
          {
          sum += neighOffset;  // sum all offsets
          }
        }
      }
    } // end of neighbor iteration

  bool pass = false;
  for( unsigned int i = 0; i < CostImageDimension; i++ )
    {
    if( ( sum[i] > 0.0 ) || ( sum[i] < 0.0 ) )
      {
      pass = true;
      }
    }

  if( !pass )
    {
    normal = newSum.normalize(); // when normal neighborhood is 0, take least
                                 // vcl_cost direction
    // std::cout << "normal is 0 " << normal << std::endl;
    }
  else
    {
    normal = sum.normalize();
    }

  /* DerivativeCostFunction = n(r)*CostValue,
      where "n(r)" is the direction of the initial point towards the lowest vcl_cost neighbor */

  derivative = DerivativeType( CostImageDimension);
  for( unsigned int i = 0; i < CostImageDimension; i++ )
    {
    derivative[i] = static_cast<double>( normal[i] * value );
    }
}
} // end of namespace itk

#endif
