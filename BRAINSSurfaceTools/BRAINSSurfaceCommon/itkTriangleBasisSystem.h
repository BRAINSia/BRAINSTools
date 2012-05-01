/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkTriangleBasisSystem.h,v $
  Language:  C++
  Date:      $Date: 2008-10-17 13:35:26 $
  Version:   $Revision: 1.47 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkTriangleBasisSystem_h
#define __itkTriangleBasisSystem_h

#include "itkMesh.h"
#include "itkPoint.h"

namespace itk
{
/** \class TriangleBasisSystem
 * \brief  Stores basis coefficients at a triangular cell.
 *
 * TriangleBasisSystem stores basis coefficients within a triangle. Basis
 * coefficients can be used thereafter for interpolation and gradient
 * computation within that triangle.
 *
 * This class is templated over the input vector type and dimension of basis.
 *
 *
 * \sa Cell
 *
 * \ingroup TriangleBasisSystems
 */
template <class TVector, unsigned int NSubspaceDimension>
class ITK_EXPORT TriangleBasisSystem
{
public:
  typedef TVector VectorType;

  /** Set/Get the vector at index k. */
  void SetVector( unsigned int k, const VectorType & v );

  const VectorType & GetVector( unsigned int k ) const;

  TriangleBasisSystem();
  virtual ~TriangleBasisSystem();

  TriangleBasisSystem( const TriangleBasisSystem & rhs );
  const TriangleBasisSystem & operator=( const TriangleBasisSystem & rhs );

private:
  VectorType m_Basis[NSubspaceDimension];
};
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkTriangleBasisSystem.txx"
#endif

#endif
