/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkTriangleBasisSystemCalculator.txx,v $
  Language:  C++
  Date:      $Date:  $
  Version:   $Revision:  $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkTriangleBasisSystemCalculator_hxx
#define __itkTriangleBasisSystemCalculator_hxx

#include "itkTriangleBasisSystemCalculator.h"
#include "itkPointLocator2.h"

namespace itk
{
/**
 * Constructor
 */
template <class TMesh, class TBasisSystem>
TriangleBasisSystemCalculator<TMesh, TBasisSystem>
::TriangleBasisSystemCalculator()
{
  itkDebugMacro("Constructor");
  this->m_InputMesh = NULL;
}

/**
 * Destructor
 */
template <class TMesh, class TBasisSystem>
TriangleBasisSystemCalculator<TMesh, TBasisSystem>
::~TriangleBasisSystemCalculator()
{
  itkDebugMacro("Destructor");
}

template <class TMesh, class TBasisSystem>
void
TriangleBasisSystemCalculator<TMesh, TBasisSystem>
::CalculateTriangle( unsigned int cellIndex, TBasisSystem & bs ) const
{
  if( this->m_InputMesh.IsNull() )
    {
    itkExceptionMacro(<< "TriangleBasisSystemCalculator CalculateTriangle  m_InputMesh is NULL.");
    }

  if( cellIndex >= this->m_InputMesh->GetNumberOfCells() )
    {
    itkExceptionMacro(<< "TriangleBasisSystemCalculator CalculateTriangle  cellIndex is too high.");
    }

  const PointsContainer * points = this->m_InputMesh->GetPoints();
  const CellsContainer *  cells = this->m_InputMesh->GetCells();

  const CellType * cell = cells->ElementAt( cellIndex );

  typedef typename CellType::PointIdentifier PointIdentifier;

  const PointIdentifier * pointIds = cell->GetPointIds();

  //
  // Get the vertexes of this triangle
  //
  PointType pt1 = points->GetElement( pointIds[0] );
  PointType pt2 = points->GetElement( pointIds[1] );
  PointType pt3 = points->GetElement( pointIds[2] );

  this->CalculateBasis( pt1, pt2, pt3, bs );
}

template <class TMesh, class TBasisSystem>
void
TriangleBasisSystemCalculator<TMesh, TBasisSystem>
::CalculateBasis(PointType pt1, PointType pt2, PointType pt3, TBasisSystem & bs ) const
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
  const double dotproduct =  v12 * v32;
  VectorType   u12 = v12 - v32 * ( dotproduct / v32.GetSquaredNorm() );
  VectorType   u32 = v32 - v12 * ( dotproduct / v12.GetSquaredNorm() );

  //
  // Add normalizations for making {u12,u32} a vector basis orthonormal to {v12, v32}.
  //
  u12 /= ( u12 * v12 );
  u32 /= ( u32 * v32 );

  bs.SetVector( 0, u12 );
  bs.SetVector( 1, u32 );
}

template <class TMesh, class TBasisSystem>
void
TriangleBasisSystemCalculator<TMesh, TBasisSystem>
::CalculateTriangle( unsigned int cellIndex, TBasisSystem & bs, TBasisSystem & bt ) const
{
  if( this->m_InputMesh.IsNull() )
    {
    itkExceptionMacro(<< "TriangleBasisSystemCalculator CalculateTriangle  m_InputMesh is NULL.");
    }

  if( cellIndex >= this->m_InputMesh->GetNumberOfCells() )
    {
    itkExceptionMacro(<< "TriangleBasisSystemCalculator CalculateTriangle  cellIndex is too high.");
    }

  const PointsContainer * points = this->m_InputMesh->GetPoints();
  const CellsContainer *  cells = this->m_InputMesh->GetCells();

  const CellType * cell = cells->ElementAt( cellIndex );

  typedef typename CellType::PointIdentifier PointIdentifier;

  const PointIdentifier * pointIds = cell->GetPointIds();

  //
  // Get the vertexes of this triangle
  //
  PointType pt1 = points->GetElement( pointIds[0] );
  PointType pt2 = points->GetElement( pointIds[1] );
  PointType pt3 = points->GetElement( pointIds[2] );

  this->CalculateBasis( pt1, pt2, pt3, bs, bt );
}

template <class TMesh, class TBasisSystem>
void
TriangleBasisSystemCalculator<TMesh, TBasisSystem>
::CalculateBasis(PointType pt1, PointType pt2, PointType pt3,
                 TBasisSystem & bs, TBasisSystem & bt ) const
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
  const double dotproduct =  v12 * v32;
  VectorType   u12 = v12 - v32 * ( dotproduct / v32.GetSquaredNorm() );
  VectorType   u32 = v32 - v12 * ( dotproduct / v12.GetSquaredNorm() );

  //
  // Add normalizations for making {u12,u32} a vector basis orthonormal to {v12, v32}.
  //
  u12 /= ( u12 * v12 );
  u32 /= ( u32 * v32 );

  bs.SetVector( 0, u12 );
  bs.SetVector( 1, u32 );

  bt.SetVector( 0, v12 );
  bt.SetVector( 1, v32 );
}
} // end namespace itk

#endif
