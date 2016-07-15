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
#ifndef __itkAssignScalarValuesQuadEdgeMeshFilter_h
#define __itkAssignScalarValuesQuadEdgeMeshFilter_h

#include "itkQuadEdgeMeshToQuadEdgeMeshFilter.h"

namespace itk
{
/**
 * \class AssignScalarValuesQuadEdgeMeshFilter
 * \brief This filter was initially proposed to assign resampled moving values back onto FixedMesh.
 *
 *
 * This filter takes two inputs as source and input meshes and assumes that
 * the points in the input are one-to-one for the points
 * in the source Mesh. Then, it assigns the scalar values of the source mesh points
 * to the points of the input mesh.
 * \ingroup MeshFilters
 *
 */
template <class TInputMesh, class TSourceMesh, class TOutputMesh>
class AssignScalarValuesQuadEdgeMeshFilter :
  public QuadEdgeMeshToQuadEdgeMeshFilter<TInputMesh, TOutputMesh>
{
public:
  typedef AssignScalarValuesQuadEdgeMeshFilter Self;
  typedef QuadEdgeMeshToQuadEdgeMeshFilter<
      TInputMesh, TOutputMesh>                     Superclass;
  typedef SmartPointer<Self>       Pointer;
  typedef SmartPointer<const Self> ConstPointer;

  /** Run-time type information (and related methods).   */
  itkTypeMacro( AssignScalarValuesQuadEdgeMeshFilter, QuadEdgeMeshToQuadEdgeMeshFilter );

  /** New macro for creation of through a Smart Pointer   */
  itkNewMacro( Self );

  typedef TSourceMesh                                  SourceMeshType;
  typedef typename  SourceMeshType::Pointer            SourceMeshPointer;
  typedef typename  SourceMeshType::PointDataContainer SourcePointDataContainer;

  typedef typename  Superclass::InputMeshType InputMeshType;

  typedef typename  Superclass::OutputMeshType           OutputMeshType;
  typedef typename  OutputMeshType::Pointer              OutputMeshPointer;
  typedef typename  Superclass::OutputPointDataContainer OutputPointDataContainer;
  typedef typename  OutputPointDataContainer::Pointer    OutputPointDataContainerPointer;
  typedef typename  OutputPointDataContainer::Iterator   OutputPointDataContainerIterator;

  // Set source mesh who has the scalar values we want.
  void SetSourceMesh( const SourceMeshType * sourceMesh );

  const SourceMeshType * GetSourceMesh( void ) const;

  // Set input mesh who wants to have scalar values from sourceMesh.
  void SetInputMesh( const InputMeshType * inputMesh );

  const InputMeshType * GetInputMesh( void ) const;

protected:
  AssignScalarValuesQuadEdgeMeshFilter();
  ~AssignScalarValuesQuadEdgeMeshFilter();

  void GenerateData() ITK_OVERRIDE;

private:

  ITK_DISALLOW_COPY_AND_ASSIGN(AssignScalarValuesQuadEdgeMeshFilter);
};
}

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkAssignScalarValuesQuadEdgeMeshFilter.hxx"
#endif

#endif
