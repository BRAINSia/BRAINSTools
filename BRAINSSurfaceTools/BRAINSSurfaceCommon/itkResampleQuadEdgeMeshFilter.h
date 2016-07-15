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
#ifndef __itkResampleQuadEdgeMeshFilter_h
#define __itkResampleQuadEdgeMeshFilter_h

#include "itkQuadEdgeMeshToQuadEdgeMeshFilter.h"
#include "itkInterpolateMeshFunction.h"
#include "itkTransform.h"

namespace itk
{
/**
 * \class ResampleQuadEdgeMeshFilter
 * \brief This resamples the scalar values of one QuadEdgeMesh into another one
 * via a user-provided Transform and Interpolator.
 *
 * \ingroup MeshFilters
 *
 */
template <class TInputMesh, class TOutputMesh>
class ResampleQuadEdgeMeshFilter :
  public QuadEdgeMeshToQuadEdgeMeshFilter<TInputMesh, TOutputMesh>
{
public:
  typedef ResampleQuadEdgeMeshFilter Self;
  typedef QuadEdgeMeshToQuadEdgeMeshFilter<
      TInputMesh, TOutputMesh>                           Superclass;
  typedef SmartPointer<Self>       Pointer;
  typedef SmartPointer<const Self> ConstPointer;

  /** Run-time type information (and related methods).   */
  itkTypeMacro( ResampleQuadEdgeMeshFilter, QuadEdgeMeshToQuadEdgeMeshFilter );

  /** New macro for creation of through a Smart Pointer   */
  itkNewMacro( Self );

  typedef TInputMesh                                 InputMeshType;
  typedef typename InputMeshType::Pointer            InputMeshPointer;
  typedef typename InputMeshType::PixelType          InputPixelType;
  typedef typename InputMeshType::PointDataContainer InputPointDataContainer;

  typedef TOutputMesh                                        OutputMeshType;
  typedef typename OutputMeshType::Pointer                   OutputMeshPointer;
  typedef typename OutputMeshType::ConstPointer              OutputMeshConstPointer;
  typedef typename OutputMeshType::EdgeCellType              OutputEdgeCellType;
  typedef typename OutputMeshType::PolygonCellType           OutputPolygonCellType;
  typedef typename OutputMeshType::PointIdList               OutputPointIdList;
  typedef typename OutputMeshType::CellTraits                OutputCellTraits;
  typedef typename OutputCellTraits::PointIdInternalIterator OutputPointsIdInternalIterator;
  typedef typename OutputMeshType::QEType                    OutputQEType;
  typedef typename OutputMeshType::PointIdentifier           OutputPointIdentifier;
  typedef typename OutputMeshType::PointType                 OutputPointType;
  typedef typename OutputPointType::VectorType               OutputVectorType;
  typedef typename OutputPointType::CoordRepType             OutputCoordType;
  typedef typename OutputMeshType::PointsContainer           OutputPointsContainer;
  typedef typename OutputMeshType::PixelType                 OutputPixelType;

  typedef typename OutputMeshType::PointsContainerConstPointer    OutputPointsContainerConstPointer;
  typedef typename OutputMeshType::PointsContainerPointer         OutputPointsContainerPointer;
  typedef typename OutputMeshType::PointsContainerIterator        OutputPointsContainerIterator;
  typedef typename OutputMeshType::PointsContainerConstIterator   OutputPointsContainerConstIterator;
  typedef typename OutputMeshType::CellsContainerPointer          OutputCellsContainerPointer;
  typedef typename OutputMeshType::CellsContainerConstPointer     OutputCellsContainerConstPointer;
  typedef typename OutputMeshType::CellsContainerIterator         OutputCellsContainerIterator;
  typedef typename OutputMeshType::CellsContainerConstIterator    OutputCellsContainerConstIterator;
  typedef typename OutputMeshType::PointDataContainer             OutputPointDataContainer;
  typedef typename OutputMeshType::PointDataContainerPointer      OutputPointDataContainerPointer;
  typedef typename OutputMeshType::PointDataContainerConstPointer OutputPointDataContainerConstPointer;
  typedef typename OutputMeshType::CellDataContainer              OutputCellDataContainer;

  itkStaticConstMacro( PointDimension, unsigned int, OutputMeshType::PointDimension );

  /** Transform typedef. */
  typedef Transform<double,
                    itkGetStaticConstMacro(PointDimension),
                    itkGetStaticConstMacro(PointDimension)>         TransformType;
  typedef typename TransformType::ConstPointer TransformPointerType;

  /** Interpolator typedef. */
  typedef InterpolateMeshFunction<InputMeshType> InterpolatorType;
  typedef typename InterpolatorType::Pointer     InterpolatorPointerType;

  /** Set Mesh whose grid will define the geometry and topology of the output Mesh.
   *  In a registration scenario, this will typically be the Fixed mesh. */
  void SetReferenceMesh( const OutputMeshType * mesh );

  const OutputMeshType * GetReferenceMesh( void ) const;

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

  /** Set the interpolator function.  The default is
   * itk::LinearInterpolateMeshFunction<InputMeshType, TInterpolatorPrecisionType>. Some
   * other options are itk::NearestNeighborInterpolateMeshFunction
   * (useful for binary masks and other images with a small number of
   * possible pixel values), and itk::BSplineInterpolateMeshFunction
   * (which provides a higher order of interpolation).  */
  itkSetObjectMacro( Interpolator, InterpolatorType );

  /** Get a pointer to the interpolator function. */
  itkGetConstObjectMacro( Interpolator, InterpolatorType );
protected:
  ResampleQuadEdgeMeshFilter();
  ~ResampleQuadEdgeMeshFilter();

  void GenerateData() ITK_OVERRIDE;

private:

  ITK_DISALLOW_COPY_AND_ASSIGN(ResampleQuadEdgeMeshFilter);

  virtual void CopyReferenceMeshToOutputMesh();

  virtual void CopyReferenceMeshToOutputMeshGeometry();

  virtual void CopyReferenceMeshToOutputMeshPoints();

  virtual void CopyReferenceMeshToOutputMeshCells();

  virtual void CopyReferenceMeshToOutputMeshEdgeCells();

  virtual void CopyReferenceMeshToOutputMeshFieldData();

  virtual void CopyReferenceMeshToOutputMeshPointData();

  virtual void CopyReferenceMeshToOutputMeshCellData();

  TransformPointerType    m_Transform;          // Coordinate transform to use
  InterpolatorPointerType m_Interpolator;       // Image function for
};
}

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkResampleQuadEdgeMeshFilter.hxx"
#endif

#endif
