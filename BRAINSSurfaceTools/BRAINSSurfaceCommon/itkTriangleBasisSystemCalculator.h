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
#ifndef __itkTriangleBasisSystemCalculator_h
#define __itkTriangleBasisSystemCalculator_h

#include "itkMesh.h"
#include "itkPoint.h"

namespace itk
{
/** \class TriangleBasisSystemCalculator
 * \brief  Computes basis coefficients at a triangular cell.
 *
 * TriangleBasisSystemCalculator computes basis coefficients
 * within a triangle. Basis coefficients can be used thereafter
 * for interpolation and gradient computation within that triangle.
 *
 * This class is templated over the input vector type and dimension of basis.
 *
 * \sa TriangleBasisSystem
 *
 */
template <class TMesh, class TBasisSystem>
class TriangleBasisSystemCalculator : public Object
{
public:
  /** Standard class typedefs. */
  typedef TriangleBasisSystemCalculator Self;
  typedef Object                        Superclass;
  typedef SmartPointer<Self>            Pointer;
  typedef SmartPointer<const Self>      ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Standard part of every itk Object. */
  itkTypeMacro(TriangleBasisSystemCalculator, Object);

  typedef TBasisSystem                       BasisSystemType;
  typedef TMesh                              MeshType;
  typedef typename MeshType::ConstPointer    MeshConstPointer;
  typedef typename MeshType::CellType        CellType;
  typedef typename MeshType::PointType       PointType;
  typedef typename MeshType::CellsContainer  CellsContainer;
  typedef typename MeshType::PointsContainer PointsContainer;
  typedef typename PointType::VectorType     VectorType;

  /** Set/Get the mesh for which the basis system will be computed. */
  itkSetConstObjectMacro( InputMesh, MeshType );
  itkGetConstObjectMacro( InputMesh, MeshType );

  /** Compute the basis system and its orthogonal at the triangular cell of the
   * Mesh that is identified by cellIndex. */
  void CalculateTriangle( unsigned int cellIndex, TBasisSystem & bs, TBasisSystem & orthogonalSystem ) const;

  /** Compute the basis system and its orthogonal given the coordinates of
   * three points. */
  void CalculateBasis(PointType pt1, PointType pt2, PointType pt3, TBasisSystem & basisSystem,
                      TBasisSystem & orthogonalSystem ) const;

  /** Compute the basis system at the triangular cell of the
   * Mesh that is identified by cellIndex. */
  void CalculateTriangle( unsigned int cellIndex, TBasisSystem & bs ) const;

  /** Compute the basis system given the coordinates of three points. */
  void CalculateBasis(PointType pt1, PointType pt2, PointType pt3, TBasisSystem & basisSystem ) const;

protected:
  TriangleBasisSystemCalculator();
  virtual ~TriangleBasisSystemCalculator();
private:
  MeshConstPointer m_InputMesh;
};
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkTriangleBasisSystemCalculator.hxx"
#endif

#endif
