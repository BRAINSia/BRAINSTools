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
#ifndef __itkTriangleBasisSystemCalculator_hxx
#define __itkTriangleBasisSystemCalculator_hxx

#include "itkTriangleBasisSystemCalculator.h"
#include "itkPointLocator2.h"

namespace itk
{
/**
 * Constructor
 */
template <typename TMesh, typename TBasisSystem>
TriangleBasisSystemCalculator<TMesh, TBasisSystem>::TriangleBasisSystemCalculator()
{
  itkDebugMacro("Constructor");
  this->m_InputMesh = nullptr;
}

/**
 * Destructor
 */
template <typename TMesh, typename TBasisSystem>
TriangleBasisSystemCalculator<TMesh, TBasisSystem>::~TriangleBasisSystemCalculator()
{
  itkDebugMacro("Destructor");
}

template <typename TMesh, typename TBasisSystem>
void
TriangleBasisSystemCalculator<TMesh, TBasisSystem>::CalculateTriangle(unsigned int cellIndex, TBasisSystem & bs) const
{
  if (this->m_InputMesh.IsNull())
  {
    itkExceptionMacro(<< "TriangleBasisSystemCalculator CalculateTriangle  m_InputMesh is NULL.");
  }

  if (cellIndex >= this->m_InputMesh->GetNumberOfCells())
  {
    itkExceptionMacro(<< "TriangleBasisSystemCalculator CalculateTriangle  cellIndex is too high.");
  }

  const PointsContainer * points = this->m_InputMesh->GetPoints();
  const CellsContainer *  cells = this->m_InputMesh->GetCells();

  const CellType * cell = cells->ElementAt(cellIndex);

  using PointIdentifier = typename CellType::PointIdentifier;

  const PointIdentifier * pointIds = cell->GetPointIds();

  //
  // Get the vertexes of this triangle
  //
  PointType pt1 = points->GetElement(pointIds[0]);
  PointType pt2 = points->GetElement(pointIds[1]);
  PointType pt3 = points->GetElement(pointIds[2]);

  this->CalculateBasis(pt1, pt2, pt3, bs);
}

template <typename TMesh, typename TBasisSystem>
void
TriangleBasisSystemCalculator<TMesh, TBasisSystem>::CalculateBasis(PointType      pt1,
                                                                   PointType      pt2,
                                                                   PointType      pt3,
                                                                   TBasisSystem & bs) const
{
  //
  // Compute Vectors along the edges.
  // These two vectors form a vector base for the 2D space of the triangle cell.
  //
  VectorType v12 = pt1 - pt2;
  VectorType v32 = pt3 - pt2;

  //
  // Compute Vectors in the dual vector base inside the 2D space of the triangle cell.
  // u12 is orthogonal to v32
  // u32 is orthogonal to v12
  //
  const double dotproduct = v12 * v32;
  VectorType   u12 = v12 - v32 * (dotproduct / v32.GetSquaredNorm());
  VectorType   u32 = v32 - v12 * (dotproduct / v12.GetSquaredNorm());

  //
  // Add normalizations for making {u12,u32} a vector basis orthonormal to {v12, v32}.
  //
  u12 /= (u12 * v12);
  u32 /= (u32 * v32);

  bs.SetVector(0, u12);
  bs.SetVector(1, u32);
}

template <typename TMesh, typename TBasisSystem>
void
TriangleBasisSystemCalculator<TMesh, TBasisSystem>::CalculateTriangle(unsigned int   cellIndex,
                                                                      TBasisSystem & bs,
                                                                      TBasisSystem & bt) const
{
  if (this->m_InputMesh.IsNull())
  {
    itkExceptionMacro(<< "TriangleBasisSystemCalculator CalculateTriangle  m_InputMesh is NULL.");
  }

  if (cellIndex >= this->m_InputMesh->GetNumberOfCells())
  {
    itkExceptionMacro(<< "TriangleBasisSystemCalculator CalculateTriangle  cellIndex is too high.");
  }

  const PointsContainer * points = this->m_InputMesh->GetPoints();
  const CellsContainer *  cells = this->m_InputMesh->GetCells();

  const CellType * cell = cells->ElementAt(cellIndex);

  using PointIdentifier = typename CellType::PointIdentifier;

  const PointIdentifier * pointIds = cell->GetPointIds();

  //
  // Get the vertexes of this triangle
  //
  PointType pt1 = points->GetElement(pointIds[0]);
  PointType pt2 = points->GetElement(pointIds[1]);
  PointType pt3 = points->GetElement(pointIds[2]);

  this->CalculateBasis(pt1, pt2, pt3, bs, bt);
}

template <typename TMesh, typename TBasisSystem>
void
TriangleBasisSystemCalculator<TMesh, TBasisSystem>::CalculateBasis(PointType      pt1,
                                                                   PointType      pt2,
                                                                   PointType      pt3,
                                                                   TBasisSystem & bs,
                                                                   TBasisSystem & bt) const
{
  //
  // Compute Vectors along the edges.
  // These two vectors form a vector base for the 2D space of the triangle cell.
  //
  VectorType v12 = pt1 - pt2;
  VectorType v32 = pt3 - pt2;

  //
  // Compute Vectors in the dual vector base inside the 2D space of the triangle cell.
  // u12 is orthogonal to v32
  // u32 is orthogonal to v12
  //
  const double dotproduct = v12 * v32;
  VectorType   u12 = v12 - v32 * (dotproduct / v32.GetSquaredNorm());
  VectorType   u32 = v32 - v12 * (dotproduct / v12.GetSquaredNorm());

  //
  // Add normalizations for making {u12,u32} a vector basis orthonormal to {v12, v32}.
  //
  u12 /= (u12 * v12);
  u32 /= (u32 * v32);

  bs.SetVector(0, u12);
  bs.SetVector(1, u32);

  bt.SetVector(0, v12);
  bt.SetVector(1, v32);
}
} // end namespace itk

#endif
