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
#ifndef __itkDeformationFieldFromTransformMeshFilter_hxx
#define __itkDeformationFieldFromTransformMeshFilter_hxx

#include "itkDeformationFieldFromTransformMeshFilter.h"
#include "itkProgressReporter.h"

namespace itk
{
template <typename TInputMesh, typename TOutputMesh>
DeformationFieldFromTransformMeshFilter<TInputMesh, TOutputMesh>
::DeformationFieldFromTransformMeshFilter()
{
  this->m_Transform = NULL;
}

template <typename TInputMesh, typename TOutputMesh>
DeformationFieldFromTransformMeshFilter<TInputMesh, TOutputMesh>
::~DeformationFieldFromTransformMeshFilter()
{
}

template <typename TInputMesh, typename TOutputMesh>
void
DeformationFieldFromTransformMeshFilter<TInputMesh, TOutputMesh>
::GenerateOutputInformation()
{
  InputMeshConstPointer inputMesh      =  this->GetInput();
  OutputMeshPointer     outputMesh     =  this->GetOutput();

  if( !inputMesh )
    {
    itkExceptionMacro(<< "Missing Input Mesh");
    }

  if( !outputMesh )
    {
    itkExceptionMacro(<< "Missing Output Mesh");
    }

  outputMesh->SetRequestedRegion( inputMesh->GetRequestedRegion() );
  outputMesh->SetBufferedRegion( inputMesh->GetBufferedRegion() );
  outputMesh->SetRequestedRegionToLargestPossibleRegion();
}

template <typename TInputMesh, typename TOutputMesh>
void
DeformationFieldFromTransformMeshFilter<TInputMesh, TOutputMesh>
::GenerateData()
{
  using InputPointsContainerPointer = typename TInputMesh::PointsContainerConstPointer;
  using OutputPointsContainerPointer = typename TOutputMesh::PointsContainerPointer;

  InputMeshConstPointer inputMesh      =  this->GetInput();
  OutputMeshPointer     outputMesh     =  this->GetOutput();

  if( !inputMesh )
    {
    itkExceptionMacro(<< "Missing Input Mesh");
    }

  if( !outputMesh )
    {
    itkExceptionMacro(<< "Missing Output Mesh");
    }

  if( this->m_Transform.IsNull() )
    {
    itkExceptionMacro(<< "Transform is not Set");
    }

  outputMesh->SetBufferedRegion( outputMesh->GetRequestedRegion() );

  InputPointsContainerPointer  inPoints  = inputMesh->GetPoints();
  OutputPointsContainerPointer outPoints = outputMesh->GetPoints();

  outPoints->Reserve( inputMesh->GetNumberOfPoints() );
  outPoints->Squeeze();  // in case the previous mesh had
                         // allocated a larger memory

  const unsigned int numberOfPoints = outputMesh->GetNumberOfPoints();

  ProgressReporter progress(this, 0, numberOfPoints);

  InputPointsContainerConstIterator inputPoint  = inPoints->Begin();
  OutputPointsContainerIterator     outputPoint = outPoints->Begin();

  while( inputPoint != inPoints->End() )
    {
    outputPoint.Value() = m_Transform->TransformPoint( inputPoint.Value() );

    progress.CompletedPixel();

    ++inputPoint;
    ++outputPoint;
    }
}

template <typename TInputMesh, typename TOutputMesh>
void
DeformationFieldFromTransformMeshFilter<TInputMesh, TOutputMesh>
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);

  if( this->m_Transform )
    {
    os << indent << "Transform: " << this->m_Transform << std::endl;
    }
}
} // end namespace itk

#endif
