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
#ifndef __itkNearestNeighborInterpolateMeshFunction_hxx
#define __itkNearestNeighborInterpolateMeshFunction_hxx

#include "itkNearestNeighborInterpolateMeshFunction.h"

namespace itk
{
/**
 * Constructor
 */
template <typename TInputMesh>
NearestNeighborInterpolateMeshFunction<TInputMesh>::NearestNeighborInterpolateMeshFunction()
{}

/**
 * Destructor
 */
template <typename TInputMesh>
NearestNeighborInterpolateMeshFunction<TInputMesh>::~NearestNeighborInterpolateMeshFunction()
{}

/**
 * Standard "PrintSelf" method
 */
template <typename TInputMesh>
void
NearestNeighborInterpolateMeshFunction<TInputMesh>::PrintSelf(std::ostream & os, Indent indent) const
{
  this->Superclass::PrintSelf(os, indent);
}

/**
 * Evaluate the mesh at a given point position.
 */
template <typename TInputMesh>
void
NearestNeighborInterpolateMeshFunction<TInputMesh>::EvaluateDerivative(const PointType & itkNotUsed(point),
                                                                       DerivativeType &  itkNotUsed(derivative)) const
{}

/**
 * Evaluate the mesh at a given point position.
 */
template <typename TInputMesh>
typename NearestNeighborInterpolateMeshFunction<TInputMesh>::OutputType
NearestNeighborInterpolateMeshFunction<TInputMesh>::Evaluate(const PointType & point) const
{
  constexpr unsigned int       numberOfNeighbors = 1;
  InstanceIdentifierVectorType result;

  this->Search(point, numberOfNeighbors, result);

  PixelType pixelValue = itk::NumericTraits<PixelType>::ZeroValue();

  const PointIdentifier pointId = result[0];

  this->GetPointData(pointId, &pixelValue);

  OutputType returnValue = static_cast<OutputType>(pixelValue);

  return returnValue;
}
} // end namespace itk

#endif
