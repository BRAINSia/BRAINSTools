/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkTriangleListBasisSystemCalculator.txx,v $
  Language:  C++
  Date:      $Date:  $
  Version:   $Revision:  $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkTriangleListBasisSystemCalculator_hxx
#define __itkTriangleListBasisSystemCalculator_hxx

#include "itkTriangleCell.h"
#include "itkTriangleBasisSystem.h"
#include "itkTriangleListBasisSystemCalculator.h"
// #include "itkPointLocator2.h"

namespace itk
{
/**
 * Constructor
 */
template <class TMesh, class TBasisSystem>
TriangleListBasisSystemCalculator<TMesh, TBasisSystem>
::TriangleListBasisSystemCalculator()
{
  itkDebugMacro("Constructor");
  this->m_InputMesh = NULL;
  this->m_BasisSystemList = BasisSystemListType::New();
}

/**
 * Destructor
 */
template <class TMesh, class TBasisSystem>
TriangleListBasisSystemCalculator<TMesh, TBasisSystem>
::~TriangleListBasisSystemCalculator()
{
  itkDebugMacro("Destructor");
}

template <class TMesh, class TBasisSystem>
void
TriangleListBasisSystemCalculator<TMesh, TBasisSystem>
::Calculate()
{
  if( this->m_InputMesh.IsNull() )
    {
    itkExceptionMacro(<< "TriangleListBasisSystemCalculator CalculateTriangle  m_InputMesh is NULL.");
    }

  this->m_BasisSystemList = BasisSystemListType::New();

  this->m_BasisSystemList->Reserve( this->m_InputMesh->GetCells()->Size() );

  typedef TriangleBasisSystemCalculator<TMesh, TBasisSystem> TriangleBasisSystemCalculatorType;

  typename TriangleBasisSystemCalculatorType::Pointer basisCalculator =
    TriangleBasisSystemCalculatorType::New();

  basisCalculator->SetInputMesh( this->m_InputMesh );

  typename CellsContainer::ConstPointer cells = this->m_InputMesh->GetCells();
  CellsContainerConstIterator cellIterator = cells->Begin();
  CellsContainerConstIterator cellEnd = cells->End();

  TBasisSystem triangleBasisSystem;

  while( cellIterator != cellEnd )
    {
    const CellIdentifier cellIndex = cellIterator.Index();
    basisCalculator->CalculateTriangle( cellIndex, this->m_BasisSystemList->ElementAt(cellIndex) );
    ++cellIterator;
    }
}
} // end namespace itk

#endif
