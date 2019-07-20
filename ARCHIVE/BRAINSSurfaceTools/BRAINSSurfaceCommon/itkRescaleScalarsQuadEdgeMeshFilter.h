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
#ifndef __itkRescaleScalarsQuadEdgeMeshFilter_h
#define __itkRescaleScalarsQuadEdgeMeshFilter_h

#include "itkQuadEdgeMeshToQuadEdgeMeshFilter.h"

namespace itk
{
/**
 * \class RescaleScalarsQuadEdgeMeshFilter
 * \brief This filter rescales the scalar values of a mesh.
 *
 * The output mesh will have the same geometry as the input.
 *
 *
 * It searches min and max values of the input mesh scalars automatically.
 * Applies a linear transformation to the scalar levels of the input Image.
 * by min and max values given by users.
 *
 * \ingroup MeshFilters
 *
 */
template < typename TMesh >
class RescaleScalarsQuadEdgeMeshFilter : public QuadEdgeMeshToQuadEdgeMeshFilter< TMesh, TMesh >
{
public:
  ITK_DISALLOW_COPY_AND_ASSIGN( RescaleScalarsQuadEdgeMeshFilter );

  using Self = RescaleScalarsQuadEdgeMeshFilter;
  using Superclass = QuadEdgeMeshToQuadEdgeMeshFilter< TMesh, TMesh >;
  using Pointer = SmartPointer< Self >;
  using ConstPointer = SmartPointer< const Self >;

  /** Run-time type information (and related methods).   */
  itkTypeMacro( RescaleScalarsQuadEdgeMeshFilter, QuadEdgeMeshToQuadEdgeMeshFilter );

  /** New macro for creation of through a Smart Pointer   */
  itkNewMacro( Self );

  using InputMeshType = TMesh;
  using InputPixelType = typename InputMeshType::PixelType;
  using InputPointDataContainer = typename InputMeshType::PointDataContainer;
  using InputPointDataContainerConstPointer = typename InputMeshType::PointDataContainerConstPointer;

  using OutputMeshType = TMesh;
  using OutputPixelType = typename OutputMeshType::PixelType;
  using OutputPointDataContainer = typename OutputMeshType::PointDataContainer;
  using OutputPointDataContainerPointer = typename OutputMeshType::PointDataContainerPointer;

  /** Set/Get the mesh that will be deformed. */
  void
  SetInputMesh( const InputMeshType * mesh );

  const InputMeshType *
  GetInputMesh( void ) const;

  /** Set/Get min value of the output mesh scalars. */
  itkSetMacro( OutputMinimum, OutputPixelType );
  itkGetMacro( OutputMinimum, OutputPixelType );

  /** Set/Get max value of the output mesh scalars. */
  itkSetMacro( OutputMaximum, OutputPixelType );
  itkGetMacro( OutputMaximum, OutputPixelType );

  /** Get min value of the input mesh scalars. */
  itkGetMacro( InputMinimum, InputPixelType );

  /** Get max value of the input mesh scalars. */
  itkGetMacro( InputMaximum, InputPixelType );

  /** Get scale to transform mesh scalars. */
  itkGetMacro( Scale, double );

protected:
  RescaleScalarsQuadEdgeMeshFilter();
  ~RescaleScalarsQuadEdgeMeshFilter();

  void
  GenerateData();

private:
  OutputPixelType m_OutputMinimum;
  OutputPixelType m_OutputMaximum;

  InputPixelType m_InputMinimum;
  InputPixelType m_InputMaximum;

  double m_Scale;
};
} // namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#  include "itkRescaleScalarsQuadEdgeMeshFilter.hxx"
#endif

#endif
