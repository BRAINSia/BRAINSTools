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
#ifndef __itkQuadEdgeMeshVectorPixelValuesSmoothingFilter_h
#define __itkQuadEdgeMeshVectorPixelValuesSmoothingFilter_h

#include "itkQuadEdgeMeshToQuadEdgeMeshFilter.h"
#include "itkQuadEdgeMeshParamMatrixCoefficients.h"

namespace itk
{
/**
 * \class QuadEdgeMeshVectorPixelValuesSmoothingFilter
 * \brief This filter smooths vector pixel values associated with points.
 *
 * This filter was based on the filter provided by
 * Arnaud Gelas, Alex Gouaillard and Sean Megason in their Insight Journal paper
 * http://hdl.handle.net/1926/1518
 * http://www.insight-journal.org/browse/publication/313
 *
 * The difference between this current filter and the one above is that this
 * filter smooths the values associated with the points (PointData) without
 * changing the actual positions of the points in space, while the filter above
 * smooths the point positions while leaving unchanged the pixel values
 * associated with the points.
 *
 * This filter expects the PixelType to be of Vector type, where the Vectors
 * will have the same number of components as the space dimension. Vectors will
 * be averaged among neighbors by first performing parallel transport to the
 * central node and then using a weighted sum. The smoothing process is performed
 * for a user-specified number of iterations.
 *
 * A full description of this filter is available in the TMI paper:
 *
 * "Spherical Demons: Fast Diffeomorphic Landmark-Free Surface Registration"
 *
 * by
 *       B.T. Thomas Yeo, Mert R. Sabuncu, Tom Vercauteren,
 *       Nicholas Ayache, Bruce Fischl, Polina Golland.
 *
 * \sa QuadEdgeMeshScalarPixelValuesSmoothingFilter
 * \ingroup MeshFilters
 *
 */
template <class TInputMesh, class TOutputMesh>
class QuadEdgeMeshVectorPixelValuesSmoothingFilter :
  public QuadEdgeMeshToQuadEdgeMeshFilter<TInputMesh, TOutputMesh>
{
public:
  typedef QuadEdgeMeshVectorPixelValuesSmoothingFilter Self;
  typedef QuadEdgeMeshToQuadEdgeMeshFilter<
      TInputMesh, TOutputMesh>                           Superclass;
  typedef SmartPointer<Self>       Pointer;
  typedef SmartPointer<const Self> ConstPointer;

  /** Run-time type information (and related methods).   */
  itkTypeMacro( QuadEdgeMeshVectorPixelValuesSmoothingFilter, QuadEdgeMeshToQuadEdgeMeshFilter );

  /** New macro for creation of through a Smart Pointer   */
  itkNewMacro( Self );

  typedef TInputMesh                                 InputMeshType;
  typedef typename InputMeshType::Pointer            InputMeshPointer;
  typedef typename InputMeshType::PixelType          InputPixelType;
  typedef typename InputMeshType::PointDataContainer InputPointDataContainer;

  typedef TOutputMesh                                        OutputMeshType;
  typedef typename OutputMeshType::Pointer                   OutputMeshPointer;
  typedef typename OutputMeshType::EdgeCellType              OutputEdgeCellType;
  typedef typename OutputMeshType::PolygonCellType           OutputPolygonCellType;
  typedef typename OutputMeshType::QEType                    OutputQEType;
  typedef typename OutputMeshType::PointIdentifier           OutputPointIdentifier;
  typedef typename OutputMeshType::PointType                 OutputPointType;
  typedef typename OutputPointType::VectorType               OutputVectorType;
  typedef typename OutputPointType::CoordRepType             OutputCoordType;
  typedef typename OutputMeshType::PointsContainer           OutputPointsContainer;
  typedef typename OutputMeshType::PointsContainerPointer    OutputPointsContainerPointer;
  typedef typename OutputMeshType::PointsContainerIterator   OutputPointsContainerIterator;
  typedef typename OutputMeshType::CellsContainerPointer     OutputCellsContainerPointer;
  typedef typename OutputMeshType::CellsContainerIterator    OutputCellsContainerIterator;
  typedef typename OutputMeshType::PointDataContainer        OutputPointDataContainer;
  typedef typename OutputMeshType::PointDataContainerPointer OutputPointDataContainerPointer;
  typedef typename OutputMeshType::PixelType                 OutputPixelType;

  itkStaticConstMacro( PointDimension, unsigned int, OutputMeshType::PointDimension );

  /** The smoothing filter will run iteratively until reaching this maximum
   * number of iterations. Emprical observartions indicate that ten iterations
   * are enough for typical deformation fields, but of course this would depend
   * on the process that you used for generating your deformation field.
   */
  itkSetMacro( MaximumNumberOfIterations, unsigned long );
  itkGetMacro( MaximumNumberOfIterations, unsigned long );

  /** Factor that controls the degree of Smoothing. Large values of Lambda
   * result is stronger smoothing.  The Lambda factor is used to compute the
   * weights of the smoothing kernel as
   *
   * \f$
   * \frac{ \exp( \frac{-1}{2 \lambda} }{ 1 + \abs{ N_i } \exp( \frac{-1}{2 \lambda} }
   * \f$
   *
   * where \f$ N_i \f$ is the number of neighbor nodes around node \f$ i \f$.
   *
   * The default value of Lambda is 1.0.
   *
   * The filter assumes that the neighbor nodes of any given nodes are located
   * at similar distances, and therefore uses the same weight for each one of
   * the neighbor values when computing their weighted average.
   *
   */
  itkSetMacro( Lambda, double );
  itkGetMacro( Lambda, double );

  /** Set Sphere Center.  The implementation of this filter assumes that the
   * Mesh surface has a spherical geometry (not only spherical topology). With
   * this method you can specify the coordinates of the center of the sphere
   * represented by the Mesh. This will be used in the computation of parallel
   * transport for vector field values associated with nodes.
   */
  itkSetMacro( SphereCenter, OutputPointType );
  itkGetConstMacro( SphereCenter, OutputPointType );

  /** Set Sphere Radius.  The implementation of this filter assumes that the
   * Mesh surface has a spherical geometry (not only spherical topology). With
   * this method you can specify the radius of the sphere. This will be used in
   * the computation of parallel transport for vector field values associated
   * with nodes.
   */
  itkSetMacro( SphereRadius, double );
  itkGetConstMacro( SphereRadius, double );
protected:
  QuadEdgeMeshVectorPixelValuesSmoothingFilter();
  ~QuadEdgeMeshVectorPixelValuesSmoothingFilter();

  void GenerateData();

  /** This method applies parallel transport to a pixel value. The typical use
   * case of this method is to manage spherical meshes whose pixel types are
   * vectors. The process of parallel transport makes possible to bring the pixel
   * values of neighbor nodes to a central node where their weighted average can
   * be computed.
   */
  void ParalelTransport(const OutputPointType src, const OutputPointType dst, const InputPixelType & inputPixelValue,
                        InputPixelType & transportedPixelValue ) const;

private:

  QuadEdgeMeshVectorPixelValuesSmoothingFilter( const Self & ); // purposely not implemented
  void operator=( const Self & );                               // purposely not implemented

  unsigned long m_MaximumNumberOfIterations;
  double        m_Lambda;

  OutputPointType m_SphereCenter;

  double m_SphereRadius;
};
}

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkQuadEdgeMeshVectorPixelValuesSmoothingFilter.txx"
#endif

#endif
