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
#ifndef __itkTriangleListBasisSystemCalculator_h
#define __itkTriangleListBasisSystemCalculator_h

#include "itkMesh.h"
#include "itkPoint.h"
#include "itkTriangleBasisSystemCalculator.h"

namespace itk
{
/** \class TriangleListBasisSystemCalculator
 * \brief  Computes basis coefficients for a list of triangular cells.
 *
 * TriangleListBasisSystemCalculator computes basis coefficients
 * within a list of triangles contained in a mesh. Basis coefficients
 * can be used thereafter for interpolation and gradient computation
 * within those triangles.
 *
 * This class is templated over the input vector type and dimension of basis.
 *
 * \sa TriangleBasisSystem
 *
 */
template < typename TMesh, typename TBasisSystem >
class TriangleListBasisSystemCalculator : public Object
{
public:
  /** Standard class type alias. */
  using Self = TriangleListBasisSystemCalculator;
  using Superclass = Object;
  using Pointer = SmartPointer< Self >;
  using ConstPointer = SmartPointer< const Self >;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Standard part of every itk Object. */
  itkTypeMacro( TriangleListBasisSystemCalculator, Object );

  using BasisSystemType = TBasisSystem;
  using CellIdentifier = typename TMesh::CellIdentifier;
  using BasisSystemListType = VectorContainer< CellIdentifier, TBasisSystem >;
  using BasisSystemListIterator = typename BasisSystemListType::ConstIterator;
  using BasisSystemListPointer = typename BasisSystemListType::Pointer;

  using MeshType = TMesh;
  using PointType = typename MeshType::PointType;
  using CellType = typename MeshType::CellType;
  using VectorType = typename PointType::VectorType;
  using MeshConstPointer = typename MeshType::ConstPointer;
  using CellsContainer = typename MeshType::CellsContainer;
  using CellsContainerIterator = typename CellsContainer::Iterator;
  using CellsContainerConstIterator = typename CellsContainer::ConstIterator;

  /** Set/Get the input mesh. */
  itkSetConstObjectMacro( InputMesh, MeshType );
  itkGetConstObjectMacro( InputMesh, MeshType );

  /** Compute the basis system at every triangle. */
  void
  Calculate();

  /** Get the list of basis systems. */
  itkGetConstObjectMacro( BasisSystemList, BasisSystemListType );

protected:
  TriangleListBasisSystemCalculator();
  virtual ~TriangleListBasisSystemCalculator();

private:
  MeshConstPointer       m_InputMesh;
  BasisSystemListPointer m_BasisSystemList;
};
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#  include "itkTriangleListBasisSystemCalculator.hxx"
#endif

#endif
