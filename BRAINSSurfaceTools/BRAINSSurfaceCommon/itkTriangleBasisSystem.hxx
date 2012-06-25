/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkTriangleBasisSystem.txx,v $
  Language:  C++
  Date:      $Date:  $
  Version:   $Revision:  $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkTriangleBasisSystem_hxx
#define __itkTriangleBasisSystem_hxx

#include "itkTriangleBasisSystem.h"

namespace itk
{
/**
 * Constructor
 */
template <class TVector, unsigned int NSubspaceDimension>
TriangleBasisSystem<TVector, NSubspaceDimension>
::TriangleBasisSystem()
{
}

/**
 * Destructor
 */
template <class TVector, unsigned int NSubspaceDimension>
TriangleBasisSystem<TVector, NSubspaceDimension>
::~TriangleBasisSystem()
{
}

/**
 * Copy Constructor
 */
template <class TVector, unsigned int NSubspaceDimension>
TriangleBasisSystem<TVector, NSubspaceDimension>
::TriangleBasisSystem( const TriangleBasisSystem & rhs )
{
  for( unsigned int i = 0; i < NSubspaceDimension; i++ )
    {
    this->m_Basis[i] = rhs.m_Basis[i];
    }
}

/**
 * Operator assignment
 */
template <class TVector, unsigned int NSubspaceDimension>
const TriangleBasisSystem<TVector, NSubspaceDimension> &
TriangleBasisSystem<TVector, NSubspaceDimension>
::operator=( const TriangleBasisSystem & rhs )
{
  for( unsigned int i = 0; i < NSubspaceDimension; i++ )
    {
    this->m_Basis[i] = rhs.m_Basis[i];
    }
  return *this;
}

template <class TVector, unsigned int NSubspaceDimension>
void
TriangleBasisSystem<TVector, NSubspaceDimension>
::SetVector( unsigned int k, const VectorType & v )
{
  if( k >= NSubspaceDimension )
    {
    itkGenericExceptionMacro(<< "TriangleBasisSystem  SetVector index k is too high.");
    }
  m_Basis[k] = v;
}

template <class TVector, unsigned int NSubspaceDimension>
const TVector &
TriangleBasisSystem<TVector, NSubspaceDimension>
::GetVector( unsigned int k ) const
{
  if( k >= NSubspaceDimension )
    {
    itkGenericExceptionMacro(<< "TriangleBasisSystem  GetVector index k is too high.");
    }
  return m_Basis[k];
}
} // end namespace itk

#endif
