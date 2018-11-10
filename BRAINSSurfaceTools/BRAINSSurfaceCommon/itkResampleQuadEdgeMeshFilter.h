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
template <typename TInputMesh, typename TOutputMesh>
class ResampleQuadEdgeMeshFilter :
  public QuadEdgeMeshToQuadEdgeMeshFilter<TInputMesh, TOutputMesh>
{
public:
  ITK_DISALLOW_COPY_AND_ASSIGN(ResampleQuadEdgeMeshFilter);

  using Self = ResampleQuadEdgeMeshFilter;
  using Superclass = QuadEdgeMeshToQuadEdgeMeshFilter<
      TInputMesh, TOutputMesh>;
  using Pointer = SmartPointer<Self>;
  using ConstPointer = SmartPointer<const Self>;

  /** Run-time type information (and related methods).   */
  itkTypeMacro( ResampleQuadEdgeMeshFilter, QuadEdgeMeshToQuadEdgeMeshFilter );

  /** New macro for creation of through a Smart Pointer   */
  itkNewMacro( Self );

  using InputMeshType = TInputMesh;
  using InputMeshPointer = typename InputMeshType::Pointer;
  using InputPixelType = typename InputMeshType::PixelType;
  using InputPointDataContainer = typename InputMeshType::PointDataContainer;

  using OutputMeshType = TOutputMesh;
  using OutputMeshPointer = typename OutputMeshType::Pointer;
  using OutputMeshConstPointer = typename OutputMeshType::ConstPointer;
  using OutputEdgeCellType = typename OutputMeshType::EdgeCellType;
  using OutputPolygonCellType = typename OutputMeshType::PolygonCellType;
  using OutputPointIdList = typename OutputMeshType::PointIdList;
  using OutputCellTraits = typename OutputMeshType::CellTraits;
  using OutputPointsIdInternalIterator = typename OutputCellTraits::PointIdInternalIterator;
  using OutputQEType = typename OutputMeshType::QEType;
  using OutputPointIdentifier = typename OutputMeshType::PointIdentifier;
  using OutputPointType = typename OutputMeshType::PointType;
  using OutputVectorType = typename OutputPointType::VectorType;
  using OutputCoordType = typename OutputPointType::CoordRepType;
  using OutputPointsContainer = typename OutputMeshType::PointsContainer;
  using OutputPixelType = typename OutputMeshType::PixelType;

  using OutputPointsContainerConstPointer = typename OutputMeshType::PointsContainerConstPointer;
  using OutputPointsContainerPointer = typename OutputMeshType::PointsContainerPointer;
  using OutputPointsContainerIterator = typename OutputMeshType::PointsContainerIterator;
  using OutputPointsContainerConstIterator = typename OutputMeshType::PointsContainerConstIterator;
  using OutputCellsContainerPointer = typename OutputMeshType::CellsContainerPointer;
  using OutputCellsContainerConstPointer = typename OutputMeshType::CellsContainerConstPointer;
  using OutputCellsContainerIterator = typename OutputMeshType::CellsContainerIterator;
  using OutputCellsContainerConstIterator = typename OutputMeshType::CellsContainerConstIterator;
  using OutputPointDataContainer = typename OutputMeshType::PointDataContainer;
  using OutputPointDataContainerPointer = typename OutputMeshType::PointDataContainerPointer;
  using OutputPointDataContainerConstPointer = typename OutputMeshType::PointDataContainerConstPointer;
  using OutputCellDataContainer = typename OutputMeshType::CellDataContainer;

  static constexpr unsigned int PointDimension = OutputMeshType::PointDimension;

  /** Transform type alias. */
  using TransformType = Transform<double,
                    Self::PointDimension,
                    Self::PointDimension>;
  using TransformPointerType = typename TransformType::ConstPointer;

  /** Interpolator type alias. */
  using InterpolatorType = InterpolateMeshFunction<InputMeshType>;
  using InterpolatorPointerType = typename InterpolatorType::Pointer;

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

  void GenerateData() override;

private:

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
