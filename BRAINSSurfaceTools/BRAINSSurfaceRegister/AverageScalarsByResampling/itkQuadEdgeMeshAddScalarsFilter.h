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
template <class TInputMesh1, class TInputMesh2, class TOutputMesh>
class QuadEdgeMeshAddScalarsFilter :
  public QuadEdgeMeshToQuadEdgeMeshFilter<TInputMesh1, TInputMesh1>
{
public:
  typedef QuadEdgeMeshAddScalarsFilter Self;
  typedef QuadEdgeMeshToQuadEdgeMeshFilter<
      TInputMesh1, TInputMesh1>                         Superclass;
  typedef SmartPointer<Self>       Pointer;
  typedef SmartPointer<const Self> ConstPointer;

  /** Run-time type information (and related methods).   */
  itkTypeMacro( QuadEdgeMeshAddScalarsFilter, QuadEdgeMeshToQuadEdgeMeshFilter );

  /** New macro for creation of through a Smart Pointer   */
  itkNewMacro( Self );

  typedef TInputMesh1 InputMeshType1;

  typedef TInputMesh2 InputMeshType2;

  typedef TOutputMesh OutputMeshType;

  /** Set/Get the input mesh 1. */
  void SetInput1( const InputMeshType1 * mesh );

  const InputMeshType1 * GetInput1( void ) const;

  /** Set/Get the input mesh2. */
  void SetInput2( const InputMeshType2 * mesh );

  const InputMeshType2 * GetInput2( void ) const;

protected:
  QuadEdgeMeshAddScalarsFilter();
  ~QuadEdgeMeshAddScalarsFilter();

  void GenerateData() ITK_OVERRIDE;

private:

  ITK_DISALLOW_COPY_AND_ASSIGN(QuadEdgeMeshAddScalarsFilter);
};
}

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkQuadEdgeMeshAddScalarsFilter.hxx"
#endif

#endif
