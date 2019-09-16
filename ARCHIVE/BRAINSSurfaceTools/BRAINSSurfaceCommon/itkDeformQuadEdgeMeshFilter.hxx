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
#ifndef __itkDeformQuadEdgeMeshFilter_hxx
#define __itkDeformQuadEdgeMeshFilter_hxx

#include "itkDeformQuadEdgeMeshFilter.h"
#include "itkProgressReporter.h"
#include "itkNumericTraitsVectorPixel.h"

namespace itk
{
template <typename TInputMesh, typename TReferenceMesh, typename TDestinationPoints>
DeformQuadEdgeMeshFilter<TInputMesh, TReferenceMesh, TDestinationPoints>::DeformQuadEdgeMeshFilter()
{
  this->SetNumberOfRequiredInputs(3);
  this->SetNumberOfRequiredOutputs(1);
  this->SetNumberOfIndexedOutputs(1);

  this->SetNthOutput(0, OutputMeshType::New());

  this->m_Interpolator = InterpolatorType::New();
  this->m_Interpolator->SetUseNearestNeighborInterpolationAsBackup(true);

  this->m_SphereRadius = 1.0;
  this->m_SphereCenter.Fill(0.0);
}

template <typename TInputMesh, typename TReferenceMesh, typename TDestinationPoints>
DeformQuadEdgeMeshFilter<TInputMesh, TReferenceMesh, TDestinationPoints>::~DeformQuadEdgeMeshFilter()
{}

template <typename TInputMesh, typename TReferenceMesh, typename TDestinationPoints>
void
DeformQuadEdgeMeshFilter<TInputMesh, TReferenceMesh, TDestinationPoints>::SetInputMesh(const InputMeshType * mesh)
{
  itkDebugMacro("setting input mesh to " << mesh);
  if (mesh != static_cast<const InputMeshType *>(this->ProcessObject::GetInput(0)))
  {
    this->ProcessObject::SetNthInput(0, const_cast<InputMeshType *>(mesh));
    this->Modified();
  }
}

template <typename TInputMesh, typename TReferenceMesh, typename TDestinationPoints>
const typename DeformQuadEdgeMeshFilter<TInputMesh, TReferenceMesh, TDestinationPoints>::InputMeshType *
DeformQuadEdgeMeshFilter<TInputMesh, TReferenceMesh, TDestinationPoints>::GetInputMesh() const
{
  Self *                surrogate = const_cast<Self *>(this);
  const InputMeshType * inputMesh = static_cast<const InputMeshType *>(surrogate->ProcessObject::GetInput(0));
  return inputMesh;
}

template <typename TInputMesh, typename TReferenceMesh, typename TDestinationPoints>
void
DeformQuadEdgeMeshFilter<TInputMesh, TReferenceMesh, TDestinationPoints>::SetReferenceMesh(
  const ReferenceMeshType * mesh)
{
  itkDebugMacro("setting input deformation mesh to " << mesh);
  if (mesh != static_cast<const ReferenceMeshType *>(this->ProcessObject::GetInput(1)))
  {
    this->ProcessObject::SetNthInput(1, const_cast<ReferenceMeshType *>(mesh));
    this->Modified();
  }
}

template <typename TInputMesh, typename TReferenceMesh, typename TDestinationPoints>
const typename DeformQuadEdgeMeshFilter<TInputMesh, TReferenceMesh, TDestinationPoints>::ReferenceMeshType *
DeformQuadEdgeMeshFilter<TInputMesh, TReferenceMesh, TDestinationPoints>::GetReferenceMesh() const
{
  Self *                    surrogate = const_cast<Self *>(this);
  const ReferenceMeshType * deformationMesh =
    static_cast<const ReferenceMeshType *>(surrogate->ProcessObject::GetInput(1));
  return deformationMesh;
}

template <typename TInputMesh, typename TReferenceMesh, typename TDestinationPoints>
void
DeformQuadEdgeMeshFilter<TInputMesh, TReferenceMesh, TDestinationPoints>::SetDestinationPoints(
  const DestinationPointsType * points)
{
  itkDebugMacro("setting input destination points to " << points);
  if (points != static_cast<const DestinationPointsType *>(this->ProcessObject::GetInput(2)))
  {
    this->ProcessObject::SetNthInput(2, const_cast<DestinationPointsType *>(points));
    this->Modified();
  }
}

template <typename TInputMesh, typename TReferenceMesh, typename TDestinationPoints>
const typename DeformQuadEdgeMeshFilter<TInputMesh, TReferenceMesh, TDestinationPoints>::DestinationPointsType *
DeformQuadEdgeMeshFilter<TInputMesh, TReferenceMesh, TDestinationPoints>::GetDestinationPoints() const
{
  Self *                        surrogate = const_cast<Self *>(this);
  const DestinationPointsType * destinationPoints =
    static_cast<const DestinationPointsType *>(surrogate->ProcessObject::GetInput(2));
  return destinationPoints;
}

template <typename TInputMesh, typename TReferenceMesh, typename TDestinationPoints>
void
DeformQuadEdgeMeshFilter<TInputMesh, TReferenceMesh, TDestinationPoints>::GenerateData()
{
  this->CopyInputMeshToOutputMesh();

  const ReferenceMeshType * referenceMesh = this->GetReferenceMesh();

  OutputMeshType * outputMesh = this->GetOutput();

  const DestinationPointsType * destinationPoints = this->GetDestinationPoints();

  const DestinationPointsContainerType * destinationPointsContainer = destinationPoints->GetPoints();

  using OutputPointsContainer = typename OutputMeshType::PointsContainer;

  OutputPointsContainer * outputPoints = outputMesh->GetPoints();

  const unsigned int numberOfPoints = outputMesh->GetNumberOfPoints();

  const unsigned int referenceNumberOfPoints = referenceMesh->GetNumberOfPoints();

  const unsigned int destinationNumberOfPoints = destinationPoints->GetNumberOfPoints();

  if (destinationNumberOfPoints != referenceNumberOfPoints)
  {
    itkExceptionMacro("Reference Mesh and Destination Points have "
                      << "different number of points " << referenceNumberOfPoints << " vs "
                      << destinationNumberOfPoints);
  }

  ProgressReporter progress(this, 0, numberOfPoints);

  this->m_Interpolator->SetSphereCenter(this->m_SphereCenter);

  this->m_Interpolator->SetInputMesh(referenceMesh);
  this->m_Interpolator->Initialize();

  using OutputPointIterator = typename OutputPointsContainer::Iterator;

  OutputPointIterator outputPointItr = outputPoints->Begin();
  OutputPointIterator outputPointEnd = outputPoints->End();

  using PointType = typename InterpolatorType::PointType;
  using VectorType = typename PointType::VectorType;

  PointType evaluatedPoint;
  PointType inputPoint;

  while (outputPointItr != outputPointEnd)
  {
    inputPoint.CastFrom(outputPointItr.Value());
    this->m_Interpolator->Evaluate(destinationPointsContainer, inputPoint, evaluatedPoint);

    //
    //  Project point to sphere surface
    //
    VectorType vectorToCenter(evaluatedPoint - this->m_SphereCenter);

    const double radialDistance = vectorToCenter.GetNorm();
    vectorToCenter *= this->m_SphereRadius / radialDistance;
    evaluatedPoint = this->m_SphereCenter + vectorToCenter;

    //
    // SetPoint() must be used here instead of a simple assignment in order to
    // maintain the topology of the QuadEdgeMesh.
    //
    outputPointItr.Value().SetPoint(evaluatedPoint);

    progress.CompletedPixel();

    ++outputPointItr;
  }
}
} // end namespace itk

#endif
