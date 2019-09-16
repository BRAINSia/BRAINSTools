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
#ifndef __itkMeshFunction_hxx
#define __itkMeshFunction_hxx

#include "itkMeshFunction.h"

namespace itk
{
/**
 * Constructor
 */
template <typename TInputMesh, typename TOutput>
MeshFunction<TInputMesh, TOutput>::MeshFunction()
{
  m_Mesh = nullptr;
}

/**
 * Standard "PrintSelf" method
 */
template <typename TInputMesh, typename TOutput>
void
MeshFunction<TInputMesh, TOutput>::PrintSelf(std::ostream & os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);

  os << indent << "InputMesh: " << m_Mesh.GetPointer() << std::endl;
}

/**
 * Initialize by setting the input mesh
 */
template <typename TInputMesh, typename TOutput>
void
MeshFunction<TInputMesh, TOutput>::SetInputMesh(const InputMeshType * ptr)
{
  this->m_Mesh = ptr;

  if (ptr)
  {
    // FIXME Add here the point locator...
  }
}

/**
 * Return the input mesh
 */
template <typename TInputMesh, typename TOutput>
const typename MeshFunction<TInputMesh, TOutput>::InputMeshType *
MeshFunction<TInputMesh, TOutput>::GetInputMesh() const
{
  return m_Mesh;
}
} // end namespace itk

#endif
