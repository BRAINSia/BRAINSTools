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
#ifndef __itkQuadEdgeMeshAddScalarsFilter_h
#define __itkQuadEdgeMeshAddScalarsFilter_h

#include "itkQuadEdgeMeshToQuadEdgeMeshFilter.h"

namespace itk
{
/**
 * \class QuadEdgeMeshAddScalarsFilter
 * \brief add scalars from two input meshes together
 *
 * This filter takes inputMesh1 and inputMesh2
 * as inputs and generates the outputMesh that has
 * scalars as the sum of that of inputMesh1 and inputMesh2.
 *
 * It requires that input meshes have the point
 * -point correspondence.
 *
 * \ingroup MeshFilters
 *
 */
template < typename TInputMesh1, typename TInputMesh2, typename TOutputMesh >
class QuadEdgeMeshAddScalarsFilter : public QuadEdgeMeshToQuadEdgeMeshFilter< TInputMesh1, TInputMesh1 >
{
public:
  ITK_DISALLOW_COPY_AND_ASSIGN( QuadEdgeMeshAddScalarsFilter );

  using Self = QuadEdgeMeshAddScalarsFilter;
  using Superclass = QuadEdgeMeshToQuadEdgeMeshFilter< TInputMesh1, TInputMesh1 >;
  using Pointer = SmartPointer< Self >;
  using ConstPointer = SmartPointer< const Self >;

  /** Run-time type information (and related methods).   */
  itkTypeMacro( QuadEdgeMeshAddScalarsFilter, QuadEdgeMeshToQuadEdgeMeshFilter );

  /** New macro for creation of through a Smart Pointer   */
  itkNewMacro( Self );

  using InputMeshType1 = TInputMesh1;

  using InputMeshType2 = TInputMesh2;

  using OutputMeshType = TOutputMesh;

  /** Set/Get the input mesh 1. */
  void
  SetInput1( const InputMeshType1 * mesh );

  const InputMeshType1 *
  GetInput1( void ) const;

  /** Set/Get the input mesh2. */
  void
  SetInput2( const InputMeshType2 * mesh );

  const InputMeshType2 *
  GetInput2( void ) const;

protected:
  QuadEdgeMeshAddScalarsFilter();
  ~QuadEdgeMeshAddScalarsFilter();

  void
  GenerateData() override;

private:
};
} // namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#  include "itkQuadEdgeMeshAddScalarsFilter.hxx"
#endif

#endif
