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
#ifndef __itkReplaceDestinationPointsQuadEdgeMeshFilter_h
#define __itkReplaceDestinationPointsQuadEdgeMeshFilter_h

#include "itkQuadEdgeMeshToQuadEdgeMeshFilter.h"

namespace itk
{
/**
 * \class ReplaceDestinationPointsQuadEdgeMeshFilter
 * \brief This filter generates replace the coordinates of the input point mesh
 * with the coordinates of the input point set.
 *
 * This filter takes as input a PointSet, and a fixed Mesh, and assumes that
 * the points in the PointSet are one-to-one destination points for the points
 * in the fixed Mesh. Then, it replaces the coordinates of input mesh points
 * with the coordinates of input point set points.
 * \ingroup MeshFilters
 *
 */
template < typename TInputMesh, typename TInputPointSet >
class ReplaceDestinationPointsQuadEdgeMeshFilter : public QuadEdgeMeshToQuadEdgeMeshFilter< TInputMesh, TInputMesh >
{
public:
  ITK_DISALLOW_COPY_AND_ASSIGN( ReplaceDestinationPointsQuadEdgeMeshFilter );

  using Self = ReplaceDestinationPointsQuadEdgeMeshFilter;
  using Superclass = QuadEdgeMeshToQuadEdgeMeshFilter< TInputMesh, TInputMesh >;
  using Pointer = SmartPointer< Self >;
  using ConstPointer = SmartPointer< const Self >;

  /** Run-time type information (and related methods).   */
  itkTypeMacro( ReplaceDestinationPointsQuadEdgeMeshFilter, QuadEdgeMeshToQuadEdgeMeshFilter );

  /** New macro for creation of through a Smart Pointer   */
  itkNewMacro( Self );

  using InputPointSetType = TInputPointSet;

  using InputMeshType = TInputMesh;

  using OutputMeshType = typename Superclass::OutputMeshType;
  using OutputPointType = typename Superclass::OutputPointType;

  /** Set Mesh whose grid defines the geometry and topology of the input PointSet.
   *  In a multi-resolution registration scenario, this will typically be the Input
   *  mesh at the current higher resolution level. */
  void
  SetInputMesh( const InputMeshType * inputMesh );

  const InputMeshType *
  GetInputMesh( void ) const;

  /** Set Mesh whose grid defines the geometry and topology of the input PointSet.
   *  In a multi-resolution registration scenario, this will typically be the Input
   *  mesh at the current higher resolution level. */
  void
  SetDestinationPoints( const InputPointSetType * destinationPointSet );

  const InputPointSetType *
  GetDestinationPoints( void ) const;

protected:
  ReplaceDestinationPointsQuadEdgeMeshFilter();
  ~ReplaceDestinationPointsQuadEdgeMeshFilter();

  void
  GenerateData();

private:
};
} // namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#  include "itkReplaceDestinationPointsQuadEdgeMeshFilter.hxx"
#endif

#endif
