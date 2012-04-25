/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkQuadEdgeMeshClampScalarsFilter.h,v $
  Language:  C++
  Date:      $Date: 2008-05-14 09:26:05 $
  Version:   $Revision: 1.21 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkQuadEdgeMeshClampScalarsFilter_txx
#define __itkQuadEdgeMeshClampScalarsFilter_txx

#include "itkQuadEdgeMeshClampScalarsFilter.h"
#include "itkProgressReporter.h"
#include "itkNumericTraitsVectorPixel.h"

namespace itk
{
template <class TInputMesh, class TOutputMesh>
QuadEdgeMeshClampScalarsFilter<TInputMesh, TOutputMesh>
::QuadEdgeMeshClampScalarsFilter()
{
  this->SetNumberOfRequiredInputs( 1 );
  this->SetNumberOfRequiredOutputs( 1 );
  this->SetNumberOfOutputs( 1 );

  this->SetNthOutput( 0, OutputMeshType::New() );

  this->m_ClampMin = false;
  this->m_ClampMax = false;
}

template <class TInputMesh, class TOutputMesh>
QuadEdgeMeshClampScalarsFilter<TInputMesh, TOutputMesh>
::~QuadEdgeMeshClampScalarsFilter()
{
}

template <class TInputMesh, class TOutputMesh>
void
QuadEdgeMeshClampScalarsFilter<TInputMesh, TOutputMesh>
::GenerateData()
{
  this->CopyInputMeshToOutputMesh();

  OutputPointDataContainerPointer outputPointData = this->GetOutput()->GetPointData();

  OutputPointDataContainerIterator outputDataItr = outputPointData->Begin();
  OutputPointDataContainerIterator outputDataEnd = outputPointData->End();

  InputPointDataContainerConstPointer inputPointData = this->GetInput()->GetPointData();

  InputPointDataContainerConstIterator inputDataItr = inputPointData->Begin();

  //
  // Clamp the output scalars when needed
  //

  outputDataItr = outputPointData->Begin();
  while( outputDataItr != outputDataEnd )
    {
    if( this->m_ClampMin || this->m_ClampMax )
      {
      // minimum end
      if( outputDataItr.Value() < this->m_OutputMinimum )
        {
        outputDataItr.Value() = this->m_OutputMinimum;
        }
      // maximum end
      if( outputDataItr.Value() > this->m_OutputMaximum )
        {
        outputDataItr.Value() = this->m_OutputMaximum;
        }
      }
    else
      {
      outputDataItr.Value() = inputDataItr.Value();
      }

    ++outputDataItr;
    ++inputDataItr;
    }
}
} // end namespace itk

#endif
