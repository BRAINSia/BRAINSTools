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
#ifndef __itkTriangleBasisSystem_hxx
#define __itkTriangleBasisSystem_hxx

#include "itkTriangleBasisSystem.h"

namespace itk
{
/**
 * Constructor
 */
template < typename TVector, unsigned int NSubspaceDimension >
TriangleBasisSystem< TVector, NSubspaceDimension >::TriangleBasisSystem()
{}

/**
 * Destructor
 */
template < typename TVector, unsigned int NSubspaceDimension >
TriangleBasisSystem< TVector, NSubspaceDimension >::~TriangleBasisSystem()
{}

/**
 * Copy Constructor
 */
template < typename TVector, unsigned int NSubspaceDimension >
TriangleBasisSystem< TVector, NSubspaceDimension >::TriangleBasisSystem( const TriangleBasisSystem & rhs )
{
  for ( unsigned int i = 0; i < NSubspaceDimension; i++ )
  {
    this->m_Basis[i] = rhs.m_Basis[i];
  }
}

/**
 * Operator assignment
 */
template < typename TVector, unsigned int NSubspaceDimension >
const TriangleBasisSystem< TVector, NSubspaceDimension > &
TriangleBasisSystem< TVector, NSubspaceDimension >::operator=( const TriangleBasisSystem & rhs )
{
  for ( unsigned int i = 0; i < NSubspaceDimension; i++ )
  {
    this->m_Basis[i] = rhs.m_Basis[i];
  }
  return *this;
}

template < typename TVector, unsigned int NSubspaceDimension >
void
TriangleBasisSystem< TVector, NSubspaceDimension >::SetVector( unsigned int k, const VectorType & v )
{
  if ( k >= NSubspaceDimension )
  {
    itkGenericExceptionMacro( << "TriangleBasisSystem  SetVector index k is too high." );
  }
  m_Basis[k] = v;
}

template < typename TVector, unsigned int NSubspaceDimension >
const TVector &
TriangleBasisSystem< TVector, NSubspaceDimension >::GetVector( unsigned int k ) const
{
  if ( k >= NSubspaceDimension )
  {
    itkGenericExceptionMacro( << "TriangleBasisSystem  GetVector index k is too high." );
  }
  return m_Basis[k];
}
} // end namespace itk

#endif
