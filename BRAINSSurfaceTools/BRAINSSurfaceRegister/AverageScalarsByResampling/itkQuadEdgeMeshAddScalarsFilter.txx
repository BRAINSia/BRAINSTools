/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkQuadEdgeMeshAddScalarsFilter.h,v $
  Language:  C++
  Date:      $Date: 2010-12-24 17:26:05 $
  Version:   $Revision: 1.21 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkQuadEdgeMeshAddScalarsFilter_txx
#define __itkQuadEdgeMeshAddScalarsFilter_txx

#include "itkQuadEdgeMeshAddScalarsFilter.h"
#include "itkProgressReporter.h"

namespace itk
{
template <class TInputMesh1, class TInputMesh2, class TOutputMesh>
QuadEdgeMeshAddScalarsFilter<TInputMesh1, TInputMesh2, TOutputMesh>
::QuadEdgeMeshAddScalarsFilter()
{
  this->SetNumberOfRequiredInputs( 2 );
  this->SetNumberOfRequiredOutputs( 1 );
  this->SetNumberOfOutputs( 1 );

  this->SetNthOutput( 0, OutputMeshType::New() );
}

template <class TInputMesh1, class TInputMesh2, class TOutputMesh>
QuadEdgeMeshAddScalarsFilter<TInputMesh1, TInputMesh2, TOutputMesh>
::~QuadEdgeMeshAddScalarsFilter()
{
}

template <class TInputMesh1, class TInputMesh2, class TOutputMesh>
void
QuadEdgeMeshAddScalarsFilter<TInputMesh1, TInputMesh2, TOutputMesh>
::SetInput1(const InputMeshType1 * mesh)
{
  itkDebugMacro("setting input1 to " << mesh);
  if( mesh != static_cast<const InputMeshType1 *>(this->ProcessObject::GetInput( 0 ) ) )
    {
    this->ProcessObject::SetNthInput(0, const_cast<InputMeshType1 *>( mesh ) );
    this->Modified();
    }
}

template <class TInputMesh1, class TInputMesh2, class TOutputMesh>
void
QuadEdgeMeshAddScalarsFilter<TInputMesh1, TInputMesh2, TOutputMesh>
::SetInput2(const InputMeshType2 * mesh)
{
  itkDebugMacro("setting input2 to " << mesh);
  if( mesh != static_cast<const InputMeshType2 *>(this->ProcessObject::GetInput( 1 ) ) )
    {
    this->ProcessObject::SetNthInput(1, const_cast<InputMeshType2 *>( mesh ) );
    this->Modified();
    }
}

template <class TInputMesh1, class TInputMesh2, class TOutputMesh>
const typename
QuadEdgeMeshAddScalarsFilter<TInputMesh1, TInputMesh2, TOutputMesh>::InputMeshType1
* QuadEdgeMeshAddScalarsFilter<TInputMesh1, TInputMesh2, TOutputMesh>
::GetInput1() const
  {
  Self *                 surrogate = const_cast<Self *>( this );
  const InputMeshType1 * inputMesh =
    static_cast<const InputMeshType1 *>( surrogate->ProcessObject::GetInput(0) );
  return inputMesh;
  }

template <class TInputMesh1, class TInputMesh2, class TOutputMesh>
const typename
QuadEdgeMeshAddScalarsFilter<TInputMesh1, TInputMesh2, TOutputMesh>::InputMeshType2
* QuadEdgeMeshAddScalarsFilter<TInputMesh1, TInputMesh2, TOutputMesh>
::GetInput2() const
  {
  Self *                 surrogate = const_cast<Self *>( this );
  const InputMeshType2 * inputMesh =
    static_cast<const InputMeshType2 *>( surrogate->ProcessObject::GetInput(1) );
  return inputMesh;
  }

template <class TInputMesh1, class TInputMesh2, class TOutputMesh>
void
QuadEdgeMeshAddScalarsFilter<TInputMesh1, TInputMesh2, TOutputMesh>
::GenerateData()
{
  this->CopyInputMeshToOutputMesh();

  const InputMeshType1 * inputMesh1 = this->GetInput1();
  const InputMeshType2 * inputMesh2 = this->GetInput2();

  OutputMeshType * outputMesh = this->GetOutput();

  const unsigned int numberOfPoints1 = inputMesh1->GetNumberOfPoints();
  const unsigned int numberOfPoints2 = inputMesh2->GetNumberOfPoints();

  if( numberOfPoints1 != numberOfPoints2 )
    {
    itkExceptionMacro("The inputs does not have the same number of points");
    }

  typedef typename InputMeshType1::PointDataContainer DisplacementVectorContainer1;
  typedef typename InputMeshType2::PointDataContainer DisplacementVectorContainer2;
  typedef typename OutputMeshType::PointDataContainer OutputDisplacementVectorContainer;

  const DisplacementVectorContainer1 * displacementVectors1 = inputMesh1->GetPointData();
  const DisplacementVectorContainer2 * displacementVectors2 = inputMesh2->GetPointData();

  OutputDisplacementVectorContainer * outputDisplacementVectors = outputMesh->GetPointData();

  typedef typename OutputDisplacementVectorContainer::Iterator OutputDisplacementVectorIterator;

  OutputDisplacementVectorIterator outputItr = outputDisplacementVectors->Begin();
  OutputDisplacementVectorIterator outputEnd = outputDisplacementVectors->End();

  typedef typename DisplacementVectorContainer1::ConstIterator InputIterator1;
  typedef typename DisplacementVectorContainer2::ConstIterator InputIterator2;

  InputIterator1 inputItr1 = displacementVectors1->Begin();
  InputIterator2 inputItr2 = displacementVectors2->Begin();

  while( outputItr != outputEnd )
    {
    outputItr.Value() = inputItr1.Value() + inputItr2.Value();

    ++outputItr;
    ++inputItr1;
    ++inputItr2;
    }
}
} // end namespace itk

#endif
