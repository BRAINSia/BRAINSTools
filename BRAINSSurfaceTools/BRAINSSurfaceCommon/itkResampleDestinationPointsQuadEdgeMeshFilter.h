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
#ifndef __itkResampleDestinationPointsQuadEdgeMeshFilter_h
#define __itkResampleDestinationPointsQuadEdgeMeshFilter_h

#include "itkMeshToMeshFilter.h"
#include "itkLinearInterpolateDeformationFieldMeshFunction.h"
#include "itkTransform.h"

namespace itk
{
/**
 * \class ResampleDestinationPointsQuadEdgeMeshFilter
 * \brief This filter resamples a collection of destination points.
 *
 * This filter takes as input a PointSet, and a fixed Mesh, and assumes that
 * the points in the PointSet are one-to-one destination points for the points
 * in the fixed Mesh. Then, it computes via linear interpolation the destination
 * points that would correspond to the locations indicated by the points of the
 * reference mesh.
 *
 * \ingroup MeshFilters
 *
 */
template <class TInputPointSet, class TFixedMesh, class TReferenceMesh, class TOutputPointSet>
class ResampleDestinationPointsQuadEdgeMeshFilter :
  public MeshToMeshFilter<TInputPointSet, TOutputPointSet>
{
public:
  typedef ResampleDestinationPointsQuadEdgeMeshFilter Self;
  typedef MeshToMeshFilter<
      TInputPointSet, TOutputPointSet>                   Superclass;
  typedef SmartPointer<Self>       Pointer;
  typedef SmartPointer<const Self> ConstPointer;

  /** Run-time type information (and related methods).   */
  itkTypeMacro( ResampleDestinationPointsQuadEdgeMeshFilter, MeshToMeshFilter );

  /** New macro for creation of through a Smart Pointer   */
  itkNewMacro( Self );

  typedef TInputPointSet                              InputPointSetType;
  typedef typename InputPointSetType::Pointer         InputPointSetPointer;
  typedef typename InputPointSetType::PointsContainer InputPointsContainer;

  typedef TFixedMesh                        FixedMeshType;
  typedef typename FixedMeshType::PointType FixedMeshPointType;

  typedef TReferenceMesh                                   ReferenceMeshType;
  typedef typename ReferenceMeshType::PointsContainer      ReferencePointsContainer;
  typedef typename ReferencePointsContainer::ConstIterator ReferencePointsContainerConstIterator;

  typedef TOutputPointSet                              OutputPointSetType;
  typedef typename OutputPointSetType::Pointer         OutputPointSetPointer;
  typedef typename OutputPointSetType::ConstPointer    OutputPointSetConstPointer;
  typedef typename OutputPointSetType::PointType       OutputPointType;
  typedef typename OutputPointSetType::PointsContainer OutputPointsContainer;

  typedef typename OutputPointSetType::PointsContainerConstPointer OutputPointsContainerConstPointer;
  typedef typename OutputPointSetType::PointsContainerPointer      OutputPointsContainerPointer;
  typedef typename OutputPointSetType::PointsContainerIterator     OutputPointsContainerIterator;

  itkStaticConstMacro( PointDimension, unsigned int, OutputPointSetType::PointDimension );

  /** Transform typedef. */
  typedef Transform<double,
                    itkGetStaticConstMacro(PointDimension),
                    itkGetStaticConstMacro(PointDimension)>         TransformType;
  typedef typename TransformType::ConstPointer TransformPointerType;

  /** Interpolator typedef. */
  typedef LinearInterpolateDeformationFieldMeshFunction<
      ReferenceMeshType, InputPointsContainer>             InterpolatorType;
  typedef typename InterpolatorType::Pointer InterpolatorPointerType;

  /** Set Mesh whose grid defines the geometry and topology of the input PointSet.
   *  In a multi-resolution registration scenario, this will typically be the Fixed
   *  mesh at the current higher resolution level. */
  void SetFixedMesh( const FixedMeshType * mesh );

  const FixedMeshType * GetFixedMesh( void ) const;

  /** Set Mesh whose grid will define the geometry of the output PointSet.
   *  In a multi-resolution registration scenario, this will typically be
   *  the Fixed mesh at the next higher resolution level. */
  void SetReferenceMesh( const ReferenceMeshType * mesh );

  const ReferenceMeshType * GetReferenceMesh( void ) const;

  /** Set the coordinate transformation.
   * Set the coordinate transform to use for resampling.  Note that this must
   * be in physical coordinates and it is the output-to-input transform, NOT
   * the input-to-output transform that you might naively expect.  By default
   * the filter uses an Identity transform. You must provide a different
   * transform here, before attempting to run the filter, if you do not want to
   * use the default Identity transform. */
  itkSetConstObjectMacro( Transform, TransformType );

  /** Get a pointer to the coordinate transform. */
  itkGetConstObjectMacro( Transform, TransformType );

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
  itkSetMacro( SphereCenter, OutputPointType );
  itkGetConstMacro( SphereCenter, OutputPointType );

  /** Set Sphere Radius.  The implementation of this class assumes that the
   * Mesh surface has a spherical geometry (not only spherical topology). With
   * this method you can specify the radius of the sphere. This will be used to
   * project destination points on the sphere after they have been interpolated.
   */
  itkSetMacro( SphereRadius, double );
  itkGetConstMacro( SphereRadius, double );
protected:
  ResampleDestinationPointsQuadEdgeMeshFilter();
  ~ResampleDestinationPointsQuadEdgeMeshFilter();

  void GenerateData() ITK_OVERRIDE;

private:

  ITK_DISALLOW_COPY_AND_ASSIGN(ResampleDestinationPointsQuadEdgeMeshFilter);

  void ProjectPointToSphereSurface( OutputPointType & point ) const;

  TransformPointerType    m_Transform;           // Coordinate transform to use
  InterpolatorPointerType m_Interpolator;        // Image function for

  OutputPointType m_SphereCenter;
  double          m_SphereRadius;
};
}

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkResampleDestinationPointsQuadEdgeMeshFilter.hxx"
#endif

#endif
