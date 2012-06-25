/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkNearestNeighborInterpolateMeshFunction.txx,v $
  Language:  C++
  Date:      $Date: 2008-10-17 13:35:26 $
  Version:   $Revision: 1.17 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkNearestNeighborInterpolateMeshFunction_hxx
#define __itkNearestNeighborInterpolateMeshFunction_hxx

#include "itkNearestNeighborInterpolateMeshFunction.h"

namespace itk
{
/**
 * Constructor
 */
template <class TInputMesh>
NearestNeighborInterpolateMeshFunction<TInputMesh>
::NearestNeighborInterpolateMeshFunction()
{
}

/**
 * Destructor
 */
template <class TInputMesh>
NearestNeighborInterpolateMeshFunction<TInputMesh>
::~NearestNeighborInterpolateMeshFunction()
{
}

/**
 * Standard "PrintSelf" method
 */
template <class TInputMesh>
void
NearestNeighborInterpolateMeshFunction<TInputMesh>
::PrintSelf( std::ostream& os, Indent indent) const
{
  this->Superclass::PrintSelf( os, indent );
}

/**
 * Evaluate the mesh at a given point position.
 */
template <class TInputMesh>
void
NearestNeighborInterpolateMeshFunction<TInputMesh>
::EvaluateDerivative( const PointType & itkNotUsed(point), DerivativeType & itkNotUsed(derivative) ) const
{
}

/**
 * Evaluate the mesh at a given point position.
 */
template <class TInputMesh>
typename
NearestNeighborInterpolateMeshFunction<TInputMesh>::OutputType
NearestNeighborInterpolateMeshFunction<TInputMesh>
::Evaluate( const PointType& point ) const
{
  const unsigned int           numberOfNeighbors = 1;
  InstanceIdentifierVectorType result;

  this->Search( point, numberOfNeighbors, result );

  PixelType pixelValue = itk::NumericTraits<PixelType>::Zero;

  const PointIdentifier pointId = result[0];

  this->GetPointData( pointId, &pixelValue );

  OutputType returnValue = static_cast<OutputType>( pixelValue );

  return returnValue;
}
} // end namespace itk

#endif
