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
#ifndef __itkDeformQuadEdgeMeshFilter_h
#define __itkDeformQuadEdgeMeshFilter_h

#include "itkQuadEdgeMeshToQuadEdgeMeshFilter.h"
#include "itkLinearInterpolateDeformationFieldMeshFunction.h"

namespace itk
{
/**
 * \class DeformQuadEdgeMeshFilter
 * \brief This filter deforms the mesh of its first input by applying the
 * deformation field implicity defined by the second and third inputs.
 *
 * This filter takes two meshes and a point set as inputs. The points of the
 * first input mesh are deformed following the deformation vectors that are
 * implied by the second mesh and the point set of destination points.  Both
 * meshes are expected to be representing a Spherical geometry with a zero
 * genus topology. Each point of the input mesh is projected onto the reference
 * mesh, and the corresponding destination points are interpolated for it.
 *
 * The user must set explicitly the values of the sphere radius and center.
 * Both meshes are expected to have the same radius and center.
 *
 * \ingroup MeshFilters
 *
 */
template <class TInputMesh, class TReferenceMesh, class TDestinationPoints>
class DeformQuadEdgeMeshFilter :
  public QuadEdgeMeshToQuadEdgeMeshFilter<TInputMesh, TInputMesh>
{
public:
  typedef DeformQuadEdgeMeshFilter                                 Self;
  typedef QuadEdgeMeshToQuadEdgeMeshFilter<TInputMesh, TInputMesh> Superclass;
  typedef SmartPointer<Self>                                       Pointer;
  typedef SmartPointer<const Self>                                 ConstPointer;

  /** Run-time type information (and related methods).   */
  itkTypeMacro( DeformQuadEdgeMeshFilter, QuadEdgeMeshToQuadEdgeMeshFilter );

  /** New macro for creation of through a Smart Pointer   */
  itkNewMacro( Self );

  typedef TInputMesh                                      InputMeshType;
  typedef TReferenceMesh                                  ReferenceMeshType;
  typedef TInputMesh                                      OutputMeshType;
  typedef TDestinationPoints                              DestinationPointsType;
  typedef typename DestinationPointsType::PointsContainer DestinationPointsContainerType;

  typedef typename InputMeshType::PointType InputPointType;

  /** Interpolator typedef. */
  typedef LinearInterpolateDeformationFieldMeshFunction<
      ReferenceMeshType, DestinationPointsContainerType>    InterpolatorType;
  typedef typename InterpolatorType::Pointer InterpolatorPointerType;

  /** Set/Get the mesh that will be deformed. */
  void SetInputMesh( const InputMeshType * mesh );

  const InputMeshType * GetInputMesh( void ) const;

  /** Set/Get the mesh that carried the deformation field as pixel data. */
  void SetReferenceMesh( const ReferenceMeshType * mesh );

  const ReferenceMeshType * GetReferenceMesh( void ) const;

  /** Set/Get the mesh that carried the deformation field as pixel data. */
  void SetDestinationPoints( const DestinationPointsType * points );

  const DestinationPointsType * GetDestinationPoints( void ) const;

  /** Set the interpolator function.  The default is a linear interpolator. */
  itkSetObjectMacro( Interpolator, InterpolatorType );

  /** Get a pointer to the interpolator function. */
  itkGetConstObjectMacro( Interpolator, InterpolatorType );

  /** Set Sphere Center.  The implementation of this class assumes that the
   * Mesh surface has a spherical geometry (not only spherical topology). With
   * this method you can specify the coordinates of the center of the sphere
   * represented by the Mesh. This will be used to project destination points
   * on the sphere after they have been interpolated.
   */
  itkSetMacro( SphereCenter, InputPointType );
  itkGetConstMacro( SphereCenter, InputPointType );

  /** Set Sphere Radius.  The implementation of this class assumes that the
   * Mesh surface has a spherical geometry (not only spherical topology). With
   * this method you can specify the radius of the sphere. This will be used to
   * project destination points on the sphere after they have been interpolated.
   */
  itkSetMacro( SphereRadius, double );
  itkGetConstMacro( SphereRadius, double );
protected:
  DeformQuadEdgeMeshFilter();
  ~DeformQuadEdgeMeshFilter();

  void GenerateData() ITK_OVERRIDE;

private:

  ITK_DISALLOW_COPY_AND_ASSIGN(DeformQuadEdgeMeshFilter);

  InterpolatorPointerType m_Interpolator;

  InputPointType m_SphereCenter;
  double         m_SphereRadius;
};
}

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkDeformQuadEdgeMeshFilter.hxx"
#endif

#endif
