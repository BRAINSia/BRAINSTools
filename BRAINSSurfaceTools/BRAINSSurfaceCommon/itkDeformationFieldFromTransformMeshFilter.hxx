/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkNodeScalarGradientCalculator.h,v $
  Language:  C++
  Date:      $Date: 2008-05-14 09:26:05 $
  Version:   $Revision: 1.21 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkDeformationFieldFromTransformMeshFilter_hxx
#define __itkDeformationFieldFromTransformMeshFilter_hxx

#include "itkDeformationFieldFromTransformMeshFilter.h"
#include "itkProgressReporter.h"

namespace itk
{
template <class TInputMesh, class TOutputMesh>
DeformationFieldFromTransformMeshFilter<TInputMesh, TOutputMesh>
::DeformationFieldFromTransformMeshFilter()
{
  this->m_Transform = NULL;
}

template <class TInputMesh, class TOutputMesh>
DeformationFieldFromTransformMeshFilter<TInputMesh, TOutputMesh>
::~DeformationFieldFromTransformMeshFilter()
{
}

template <class TInputMesh, class TOutputMesh>
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

template <class TInputMesh, class TOutputMesh>
void
DeformationFieldFromTransformMeshFilter<TInputMesh, TOutputMesh>
::GenerateData()
{
  typedef typename TInputMesh::PointsContainer  InputPointsContainer;
  typedef typename TOutputMesh::PointsContainer OutputPointsContainer;

  typedef typename TInputMesh::PointsContainerConstPointer InputPointsContainerPointer;
  typedef typename TOutputMesh::PointsContainerPointer     OutputPointsContainerPointer;

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

template <class TInputMesh, class TOutputMesh>
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
