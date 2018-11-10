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
template <typename TVector, unsigned int NSubspaceDimension>
class TriangleBasisSystem
{
public:
  using VectorType = TVector;

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
#include "itkTriangleBasisSystem.hxx"
#endif

#endif
