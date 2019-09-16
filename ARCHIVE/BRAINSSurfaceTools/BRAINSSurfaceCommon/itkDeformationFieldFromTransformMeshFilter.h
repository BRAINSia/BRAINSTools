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
#ifndef __itkDeformationFieldFromTransformMeshFilter_h
#define __itkDeformationFieldFromTransformMeshFilter_h

#include "itkMeshToMeshFilter.h"
#include "itkTransform.h"

namespace itk
{
/**
 * \class DeformationFieldFromTransformMeshFilter
 * \brief Generate destination points from a Mesh and a Transform.
 *
 * This filter takes as input a Mesh and a Transform and produces as
 * output a point container with the list of Mesh points image resulting from
 * the Transform mapping.
 *
 * \ingroup MeshFilters
 *
 */
template <typename TInputMesh, typename TOutputMesh>
class DeformationFieldFromTransformMeshFilter : public MeshToMeshFilter<TInputMesh, TOutputMesh>
{
public:
  ITK_DISALLOW_COPY_AND_ASSIGN(DeformationFieldFromTransformMeshFilter);

  using Self = DeformationFieldFromTransformMeshFilter;
  using Superclass = MeshToMeshFilter<TInputMesh, TOutputMesh>;
  using Pointer = SmartPointer<Self>;
  using ConstPointer = SmartPointer<const Self>;

  /** Run-time type information (and related methods).   */
  itkTypeMacro(DeformationFieldFromTransformMeshFilter, MeshToMeshFilter);

  /** New macro for creation of through a Smart Pointer   */
  itkNewMacro(Self);

  using InputMeshType = TInputMesh;
  using InputMeshConstPointer = typename InputMeshType::ConstPointer;
  using InputPointsContainerConstIterator = typename InputMeshType::PointsContainerConstIterator;

  using OutputMeshType = TOutputMesh;
  using OutputMeshPointer = typename OutputMeshType::Pointer;

  using OutputPointsContainerIterator = typename OutputMeshType::PointsContainerIterator;
  using OutputPointsContainerConstIterator = typename OutputMeshType::PointsContainerConstIterator;

  static constexpr unsigned int PointDimension = OutputMeshType::PointDimension;

  /** Transform type alias. */
  using TransformType = Transform<double, Self::PointDimension, Self::PointDimension>;
  using TransformPointerType = typename TransformType::ConstPointer;

  /** Set the coordinate transformation.  Set the coordinate transform that
   * will map the points of the input mesh to points of the output PointSet.
   * The points of the output PointSet are one-to-one the result of taking
   * points from the input Mesh and mapping them through the Transform.
   */
  itkSetConstObjectMacro(Transform, TransformType);

  /** Get a pointer to the coordinate transform. */
  itkGetConstObjectMacro(Transform, TransformType);

protected:
  DeformationFieldFromTransformMeshFilter();
  ~DeformationFieldFromTransformMeshFilter();

  void
  PrintSelf(std::ostream & os, Indent indent) const;

  void
  GenerateOutputInformation();

  void
  GenerateData();

private:
  TransformPointerType m_Transform; // Coordinate transform to use
};
} // namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#  include "itkDeformationFieldFromTransformMeshFilter.hxx"
#endif

#endif
