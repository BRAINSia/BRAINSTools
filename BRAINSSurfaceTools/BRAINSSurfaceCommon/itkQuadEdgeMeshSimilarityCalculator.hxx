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
#ifndef __itkQuadEdgeMeshSimilarityCalculator_hxx
#define __itkQuadEdgeMeshSimilarityCalculator_hxx

#include "itkQuadEdgeMeshSimilarityCalculator.h"
#include "itkProgressReporter.h"

namespace itk
{
template <class TInputMesh1, class TInputMesh2>
QuadEdgeMeshSimilarityCalculator<TInputMesh1, TInputMesh2>
::QuadEdgeMeshSimilarityCalculator()
{
  this->SetNumberOfRequiredInputs( 2 );
  this->SetNumberOfRequiredOutputs( 0 );
  this->SetNumberOfIndexedOutputs( 0 );

  this->m_Dice = 0.0;
  this->m_Jaccard = 0.0;
}

template <class TInputMesh1, class TInputMesh2>
QuadEdgeMeshSimilarityCalculator<TInputMesh1, TInputMesh2>
::~QuadEdgeMeshSimilarityCalculator()
{
}

template <class TInputMesh1, class TInputMesh2>
void
QuadEdgeMeshSimilarityCalculator<TInputMesh1, TInputMesh2>
::SetInputMesh1( const InputMeshType1 * mesh1 )
{
  itkDebugMacro("setting input mesh 1 to " << mesh1);
  if( mesh1 != static_cast<const InputMeshType1 *>(this->ProcessObject::GetInput( 0 ) ) )
    {
    this->ProcessObject::SetNthInput(0, const_cast<InputMeshType1 *>( mesh1 ) );
    this->Modified();
    }
}

template <class TInputMesh1, class TInputMesh2>
const typename
QuadEdgeMeshSimilarityCalculator<TInputMesh1, TInputMesh2>::InputMeshType1
* QuadEdgeMeshSimilarityCalculator<TInputMesh1, TInputMesh2>
::GetInputMesh1() const
  {
  Self *                 surrogate = const_cast<Self *>( this );
  const InputMeshType1 * inputMesh1 =
    static_cast<const InputMeshType1 *>( surrogate->ProcessObject::GetInput(0) );
  return inputMesh1;
  }

template <class TInputMesh1, class TInputMesh2>
void
QuadEdgeMeshSimilarityCalculator<TInputMesh1, TInputMesh2>
::SetInputMesh2( const InputMeshType2 * mesh2 )
{
  itkDebugMacro("setting input mesh 2 to " << mesh2);
  if( mesh2 != static_cast<const InputMeshType2 *>(this->ProcessObject::GetInput(1) ) )
    {
    this->ProcessObject::SetNthInput(1, const_cast<InputMeshType2 *>( mesh2 ) );
    this->Modified();
    }
}

template <class TInputMesh1, class TInputMesh2>
const typename
QuadEdgeMeshSimilarityCalculator<TInputMesh1, TInputMesh2>::InputMeshType2
* QuadEdgeMeshSimilarityCalculator<TInputMesh1, TInputMesh2>
::GetInputMesh2() const
  {
  Self *                 surrogate = const_cast<Self *>( this );
  const InputMeshType2 * inputMesh2 =
    static_cast<const InputMeshType2 *>( surrogate->ProcessObject::GetInput(1) );
  return inputMesh2;
  }

template <class TInputMesh1, class TInputMesh2>
void
QuadEdgeMeshSimilarityCalculator<TInputMesh1, TInputMesh2>
::Compute()
{
  const InputMeshType1 * inputMesh1 = this->GetInputMesh1();
  const InputMeshType2 * inputMesh2 = this->GetInputMesh2();

  const CellsContainer1 * cells1 =  inputMesh1->GetCells();
  const CellsContainer2 * cells2 =  inputMesh2->GetCells();

  const InputPointDataContainer1 * pointDataContainer1 =  inputMesh1->GetPointData();
  const InputPointDataContainer2 * pointDataContainer2 =  inputMesh2->GetPointData();

  const unsigned int numberOfPoints1 = inputMesh1->GetNumberOfPoints();
  const unsigned int numberOfPoints2 = inputMesh2->GetNumberOfPoints();

  const unsigned int numberOfCells1 = inputMesh1->GetNumberOfCells();
  const unsigned int numberOfCells2 = inputMesh2->GetNumberOfCells();

  if( (numberOfPoints1 != numberOfPoints2) || (numberOfCells1 != numberOfCells2) )
    {
    itkExceptionMacro("InputMesh1 and InputMesh2 have "
                      << "different number of points or cells" );
    }

  const unsigned int numberOfVerticesInTriangle = 3;

  InputPointType1 point1[numberOfVerticesInTriangle];

  CellsContainerConstIterator1 cellIterator1 = cells1->Begin();
  CellsContainerConstIterator1 cellEnd1 = cells1->End();

  CellsContainerConstIterator2 cellIterator2 = cells2->Begin();

  double Dice_numerator = 0.0;
  double Dice_denominator = 0.0;

  double Jaccard_numerator = 0.0;
  double Jaccard_denominator = 0.0;

  while( cellIterator1 != cellEnd1 )
    {
    CellType1 * cellPointer1 = cellIterator1.Value();
    CellType2 * cellPointer2 = cellIterator2.Value();

    // Consider current cell. Iterate through its points.
    PointIdIterator1 pointIdIterator1 = cellPointer1->PointIdsBegin();
    PointIdIterator1 pointIdEnd1 = cellPointer1->PointIdsEnd();

    // calculate the cell's area
    unsigned int i = 0;
    while( pointIdIterator1 != pointIdEnd1 )
      {
      const PointIdentifier1 pointId = *pointIdIterator1;
      point1[i] = inputMesh1->GetPoint( pointId );
      i++;
      ++pointIdIterator1;
      }

    const AreaType area = TriangleType::ComputeArea( point1[0], point1[1], point1[2] );
    // std::cout<<area<<std::endl;

    // go through points of the cell
    pointIdIterator1 = cellPointer1->PointIdsBegin();
    PointIdIterator2 pointIdIterator2 = cellPointer2->PointIdsBegin();

    InputPixelType1 pointLabel1;
    InputPixelType2 pointLabel2;
    while( pointIdIterator1 != pointIdEnd1 )
      {
      const PointIdentifier1 pointId1 = *pointIdIterator1;
      const PointIdentifier2 pointId2 = *pointIdIterator2;

      pointLabel1 = pointDataContainer1->GetElement( pointId1 );
      pointLabel2 = pointDataContainer2->GetElement( pointId2 );

      if( pointLabel1 == this->m_LabelValue )
        {
        Dice_denominator += area * (1.0 / 3.0);
        Jaccard_denominator += area * (1.0 / 3.0);
        }

      if( pointLabel2 == static_cast<InputPixelType2>(this->m_LabelValue) )
        {
        Dice_denominator += area * (1.0 / 3.0);
        Jaccard_denominator += area * (1.0 / 3.0);
        }

      if( ( pointLabel1 == this->m_LabelValue)
          && ( pointLabel2 == static_cast<InputPixelType2>(this->m_LabelValue) ) )
        {
        Dice_numerator += area * (1.0 / 3.0);
        Jaccard_numerator += area * (1.0 / 3.0);
        Jaccard_denominator -= area * (1.0 / 3.0);
        }

      ++pointIdIterator1;
      ++pointIdIterator2;
      }

    ++cellIterator1;
    ++cellIterator2;
    } // cell

  if( Dice_denominator != 0.0 )
    {
    this->m_Dice = 2.0 * Dice_numerator / Dice_denominator;
    }
  else
    {
    itkExceptionMacro("Denominator of Dice is zero!" );
    }

  if( Jaccard_denominator != 0.0 )
    {
    this->m_Jaccard = Jaccard_numerator / Jaccard_denominator;
    }
  else
    {
    itkExceptionMacro("Denominator of Jaccard is zero!" );
    }
} // GenerateData
} // end namespace itk

#endif
