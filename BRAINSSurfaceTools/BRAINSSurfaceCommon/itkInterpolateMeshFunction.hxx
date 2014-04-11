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
