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
#ifndef __itkMeshToMeshMetric_hxx
#define __itkMeshToMeshMetric_hxx

#include "itkMeshToMeshMetric.h"

namespace itk
{
/** Constructor */
template < typename TFixedMesh, typename TMovingMesh >
MeshToMeshMetric< TFixedMesh, TMovingMesh >::MeshToMeshMetric()
{
  m_FixedMesh = nullptr;    // has to be provided by the user.
  m_MovingMesh = nullptr;   // has to be provided by the user.
  m_Transform = nullptr;    // has to be provided by the user.
  m_Interpolator = nullptr; // has to be provided by the user.
}

/** Set the parameters that define a unique transform */
template < typename TFixedMesh, typename TMovingMesh >
void
MeshToMeshMetric< TFixedMesh, TMovingMesh >::SetTransformParameters( const ParametersType & parameters ) const
{
  if ( !m_Transform )
  {
    itkExceptionMacro( << "Transform has not been assigned" );
  }
  m_Transform->SetParameters( parameters );
}

/** Initialize the metric */
template < typename TFixedMesh, typename TMovingMesh >
void
MeshToMeshMetric< TFixedMesh, TMovingMesh >::Initialize( void ) throw( ExceptionObject )
{
  if ( !m_Transform )
  {
    itkExceptionMacro( << "Transform is not present" );
  }

  if ( !m_Interpolator )
  {
    itkExceptionMacro( << "Interpolator is not present" );
  }

  if ( !m_MovingMesh )
  {
    itkExceptionMacro( << "MovingMesh is not present" );
  }

  if ( !( m_MovingMesh->GetPoints() ) )
  {
    itkExceptionMacro( << "MovingMesh does not have points" );
  }

  if ( !( m_MovingMesh->GetPointData() ) )
  {
    itkExceptionMacro( << "MovingMesh does not have point data" );
  }

  if ( !m_FixedMesh )
  {
    itkExceptionMacro( << "FixedMesh is not present" );
  }

  if ( !( m_FixedMesh->GetPoints() ) )
  {
    itkExceptionMacro( << "FixedMesh does not have points" );
  }

  if ( !( m_FixedMesh->GetPointData() ) )
  {
    itkExceptionMacro( << "FixedMesh does not have point data" );
  }

  // If the Mesh is provided by a source, update the source.
  if ( m_MovingMesh->GetSource() )
  {
    m_MovingMesh->GetSource()->Update();
  }

  // If the point set is provided by a source, update the source.
  if ( m_FixedMesh->GetSource() )
  {
    m_FixedMesh->GetSource()->Update();
  }

  // Initialize the internal point locator structure
  this->m_Interpolator->SetInputMesh( this->m_MovingMesh );
  this->m_Interpolator->Initialize();
}

/** PrintSelf */
template < typename TFixedMesh, typename TMovingMesh >
void
MeshToMeshMetric< TFixedMesh, TMovingMesh >::PrintSelf( std::ostream & os, Indent indent ) const
{
  Superclass::PrintSelf( os, indent );

  os << indent << "Moving Mesh: " << m_MovingMesh.GetPointer() << std::endl;
  os << indent << "Fixed  Mesh: " << m_FixedMesh.GetPointer() << std::endl;
  os << indent << "Transform:    " << m_Transform.GetPointer() << std::endl;
  os << indent << "Interpolator: " << m_Interpolator.GetPointer() << std::endl;
}
} // end namespace itk

#endif
