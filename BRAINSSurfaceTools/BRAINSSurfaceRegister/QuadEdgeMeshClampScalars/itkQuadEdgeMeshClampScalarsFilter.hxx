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
  this->SetNumberOfIndexedOutputs( 1 );

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

  OutputPointDataContainerIterator outputDataEnd = outputPointData->End();

  InputPointDataContainerConstPointer inputPointData = this->GetInput()->GetPointData();

  InputPointDataContainerConstIterator inputDataItr = inputPointData->Begin();

  //
  // Clamp the output scalars when needed
  //

  OutputPointDataContainerIterator outputDataItr = outputPointData->Begin();
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
