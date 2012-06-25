/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkNodeScalarGradientCalculator.h,v $
  Language:  C++
  Date:      $Date: 2008-05-14 09:26:05 $
  Version:   $Revision: 1.21 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkQuadEdgeMeshGenerateDeformationFieldFilter_h
#define __itkQuadEdgeMeshGenerateDeformationFieldFilter_h

#include "itkQuadEdgeMeshToQuadEdgeMeshFilter.h"

namespace itk
{
/**
 * \class QuadEdgeMeshGenerateDeformationFieldFilter
 * \brief This filter generates a deformation field from an input mesh and a
 * list of destination points.
 *
 * This filter takes as input a PointSet, and a fixed Mesh, and assumes that
 * the points in the PointSet are one-to-one destination points for the points
 * in the fixed Mesh. Then, it computes a vector field where every node gets
 * assigned a (tangent) vector that points to its destination point. The filter expects
 * a one-to-one correspondence between the nodes of the mesh and the points
 * in the collection of destination points.
 *
 * \ingroup MeshFilters
 *
 */
template <class TInputMesh, class TInputPointSet, class TOutputMesh>
class QuadEdgeMeshGenerateDeformationFieldFilter :
  public QuadEdgeMeshToQuadEdgeMeshFilter<TInputMesh, TOutputMesh>
{
public:
  typedef QuadEdgeMeshGenerateDeformationFieldFilter Self;
  typedef QuadEdgeMeshToQuadEdgeMeshFilter<
      TInputMesh, TOutputMesh>                           Superclass;
  typedef SmartPointer<Self>       Pointer;
  typedef SmartPointer<const Self> ConstPointer;

  /** Run-time type information (and related methods).   */
  itkTypeMacro( QuadEdgeMeshGenerateDeformationFieldFilter, QuadEdgeMeshToQuadEdgeMeshFilter );

  /** New macro for creation of through a Smart Pointer   */
  itkNewMacro( Self );

  typedef TInputPointSet InputPointSetType;

  typedef TInputMesh                        InputMeshType;
  typedef typename InputMeshType::PointType PointType;

  typedef TOutputMesh OutputMeshType;

  /** Set Mesh whose grid defines the geometry and topology of the input PointSet.
   *  In a multi-resolution registration scenario, this will typically be the Input
   *  mesh at the current higher resolution level. */
  void SetInputMesh( const InputMeshType * inputMesh );

  const InputMeshType * GetInputMesh( void ) const;

  /** Set Mesh whose grid defines the geometry and topology of the input PointSet.
   *  In a multi-resolution registration scenario, this will typically be the Input
   *  mesh at the current higher resolution level. */
  void SetDestinationPoints( const InputPointSetType * destinationPointSet );

  const InputPointSetType * GetDestinationPoints( void ) const;

  /** Set Sphere Center.  The implementation of this filter assumes that the
   * Mesh surface has a spherical geometry (not only spherical topology). With
   * this method you can specify the coordinates of the center of the sphere
   * represented by the Mesh. This will be used in the computation of parallel
   * transport for vector values associated with nodes.
   */
  itkSetMacro( SphereCenter, PointType );
  itkGetConstMacro( SphereCenter, PointType );

  /** Set Sphere Radius.  The implementation of this filter assumes that the
   * Mesh surface has a spherical geometry (not only spherical topology). With
   * this method you can specify the radius of the sphere. This will be used in
   * the computation of parallel transport for vector values associated
   * with nodes.
   */
  itkSetMacro( SphereRadius, double );
  itkGetConstMacro( SphereRadius, double );
protected:
  QuadEdgeMeshGenerateDeformationFieldFilter();
  ~QuadEdgeMeshGenerateDeformationFieldFilter();

  void GenerateData();

private:

  QuadEdgeMeshGenerateDeformationFieldFilter( const Self & ); // purposely not implemented
  void operator=( const Self & );                             // purposely not implemented

  /** Center of spherical mesh. We assume that both the Fixed and
   * Moving meshes have spherical geometry and that they share the same
   * center and radius. */
  PointType m_SphereCenter;

  /** Radius of spherical mesh. We assume that both the Fixed and
   * Moving meshes have spherical geometry and that they share the same
   * center and radius. */
  double m_SphereRadius;
};
}

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkQuadEdgeMeshGenerateDeformationFieldFilter.hxx"
#endif

#endif
