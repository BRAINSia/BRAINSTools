/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkLinearInterpolateDeformationFieldMeshFunction.txx,v $
  Language:  C++
  Date:      $Date: 2008-10-17 13:35:26 $
  Version:   $Revision: 1.17 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkLinearInterpolateDeformationFieldMeshFunction_hxx
#define __itkLinearInterpolateDeformationFieldMeshFunction_hxx

#include "itkVector.h"
#include "itkQuadEdgeMesh.h"
#include "itkLinearInterpolateDeformationFieldMeshFunction.h"
#include "itkTriangleCell.h"

namespace itk
{
/**
 * Constructor
 */
template <class TInputMesh, class TDestinationPointsContainer>
LinearInterpolateDeformationFieldMeshFunction<TInputMesh, TDestinationPointsContainer>
::LinearInterpolateDeformationFieldMeshFunction()
{
}

/**
 * Destructor
 */
template <class TInputMesh, class TDestinationPointsContainer>
LinearInterpolateDeformationFieldMeshFunction<TInputMesh, TDestinationPointsContainer>
::~LinearInterpolateDeformationFieldMeshFunction()
{
}

/**
 * Standard "PrintSelf" method
 */
template <class TInputMesh, class TDestinationPointsContainer>
void
LinearInterpolateDeformationFieldMeshFunction<TInputMesh, TDestinationPointsContainer>
::PrintSelf( std::ostream& os, Indent indent) const
{
  this->Superclass::PrintSelf( os, indent );
}

/**
 * Method provided for completness of the base class API.
 * This method is not expected to be used here.
 */
template <class TInputMesh, class TDestinationPointsContainer>
typename
LinearInterpolateDeformationFieldMeshFunction<TInputMesh, TDestinationPointsContainer>::OutputType
LinearInterpolateDeformationFieldMeshFunction<TInputMesh, TDestinationPointsContainer>
::Evaluate( const PointType & itkNotUsed(point) ) const
{
  return OutputType();
}

/**
 * Evaluate the mesh at a given point position.
 */
template <class TInputMesh, class TDestinationPointsContainer>
bool
LinearInterpolateDeformationFieldMeshFunction<TInputMesh, TDestinationPointsContainer>
::Evaluate( const DestinationPointsContainerType * field,
            const PointType & point, PointType & outputPoint ) const
{
  InstanceIdentifierVectorType pointIds(3);

  bool foundTriangle = this->FindTriangle( point, pointIds );

  if( !foundTriangle )
    {
    if( this->GetUseNearestNeighborInterpolationAsBackup() )
      {
      const unsigned int numberOfNeighbors = 1;

      InstanceIdentifierVectorType closestPointIds(numberOfNeighbors);

      this->Search( point, numberOfNeighbors, closestPointIds );

      outputPoint = field->ElementAt( closestPointIds[0] );

      return true;
      }
    else
      {
      std::cout << "no triangle is found in DF interpolation." << std::endl;
      return false;
      }
    }

  const PointType & point1 = field->ElementAt( pointIds[0] );
  const PointType & point2 = field->ElementAt( pointIds[1] );
  const PointType & point3 = field->ElementAt( pointIds[2] );

  const RealType & weight1 = this->GetInterpolationWeight(0);
  const RealType & weight2 = this->GetInterpolationWeight(1);

  outputPoint.SetToBarycentricCombination( point1, point2, point3, weight1, weight2 );

  return true;
}
} // end namespace itk

#endif
