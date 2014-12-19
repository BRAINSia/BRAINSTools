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
  this->m_InputMesh = ITK_NULLPTR;
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
