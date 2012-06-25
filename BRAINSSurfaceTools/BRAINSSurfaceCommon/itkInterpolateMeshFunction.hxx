/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkInterpolateMeshFunction.txx,v $
  Language:  C++
  Date:      $Date: 2008-10-17 13:35:26 $
  Version:   $Revision: 1.17 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkInterpolateMeshFunction_hxx
#define __itkInterpolateMeshFunction_hxx

#include "itkInterpolateMeshFunction.h"

namespace itk
{
/**
 * Constructor
 */
template <class TInputMesh>
InterpolateMeshFunction<TInputMesh>
::InterpolateMeshFunction()
{
  this->m_PointLocator = PointLocatorType::New();
}

/**
 * Destructor
 */
template <class TInputMesh>
InterpolateMeshFunction<TInputMesh>
::~InterpolateMeshFunction()
{
}

/**
 * Prepare the internal data structures of the point locator
 */
template <class TInputMesh>
void
InterpolateMeshFunction<TInputMesh>
::Initialize()
{
  this->m_PointLocator->SetPointSet( this->m_Mesh );
  this->m_PointLocator->Initialize();
}

template <class TInputMesh>
void
InterpolateMeshFunction<TInputMesh>
::Search(const PointType & query,
         unsigned int numberOfNeighborsRequested,
         InstanceIdentifierVectorType& result) const
{
  typename PointLocatorType::PointType point( query );
  this->m_PointLocator->Search( point, numberOfNeighborsRequested, result );
}

template <class TInputMesh>
void
InterpolateMeshFunction<TInputMesh>
::Search(const PointType & query,
         double radius,
         InstanceIdentifierVectorType& result) const
{
  typename PointLocatorType::PointType point( query );
  this->m_PointLocator->Search( point, radius, result );
}

/**
 * Return the pixel value by delegating to the mesh.
 */
template <class TInputMesh>
void
InterpolateMeshFunction<TInputMesh>
::GetPointData( PointIdentifier pointId, PixelType * value ) const
{
  this->m_Mesh->GetPointData( pointId, value );
}

/**
 * Standard "PrintSelf" method
 */
template <class TInputMesh>
void
InterpolateMeshFunction<TInputMesh>
::PrintSelf( std::ostream& os, Indent indent) const
{
  this->Superclass::PrintSelf( os, indent );
}
} // end namespace itk

#endif
