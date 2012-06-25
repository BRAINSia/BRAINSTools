/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkMeshToMeshMetric.h,v $
  Language:  C++
  Date:      $Date: 2003-11-08 17:58:32 $
  Version:   $Revision: 1.1 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkMeshToMeshMetric_h
#define __itkMeshToMeshMetric_h

#include "itkTransform.h"
#include "itkSingleValuedCostFunction.h"
#include "itkInterpolateMeshFunction.h"
#include "itkExceptionObject.h"
#include "itkSpatialObject.h"

namespace itk
{
/** \class MeshToMeshMetric
 * \brief Computes similarity between two meshes.
 *
 * This Class is templated over the type of the two meshes.  It expects a
 * Transform to be plugged in.  This particular class is the base class for a
 * hierarchy of mesh to mesh metrics.
 *
 * This class computes a value that measures the similarity between the fixed
 * mesh and the transformed moving mesh.
 *
 * \ingroup RegistrationMetrics
 *
 */

template <class TFixedMesh,  class TMovingMesh>
class ITK_EXPORT MeshToMeshMetric : public SingleValuedCostFunction
{
public:

  /** Standard class typedefs. */
  typedef MeshToMeshMetric         Self;
  typedef SingleValuedCostFunction Superclass;
  typedef SmartPointer<Self>       Pointer;
  typedef SmartPointer<const Self> ConstPointer;

  /** Type used for representing point components  */
  typedef typename TFixedMesh::CoordRepType CoordinateRepresentationType;

  /** Run-time type information (and related methods). */
  itkTypeMacro(MeshToMeshMetric, SingleValuedCostFunction);

  /**  Type of the moving Mesh. */
  typedef TMovingMesh                           MovingMeshType;
  typedef typename TMovingMesh::PixelType       MovingMeshPixelType;
  typedef typename MovingMeshType::ConstPointer MovingMeshConstPointer;

  /**  Type of the fixed Mesh. */
  typedef TFixedMesh                           FixedMeshType;
  typedef typename FixedMeshType::ConstPointer FixedMeshConstPointer;

  /** Constants for the pointset dimensions */
  itkStaticConstMacro(MovingMeshDimension, unsigned int,
                      TMovingMesh::PointDimension);
  itkStaticConstMacro(FixedMeshDimension, unsigned int,
                      TFixedMesh::PointDimension);

  typedef typename FixedMeshType::PointsContainer::ConstIterator    PointIterator;
  typedef typename FixedMeshType::PointDataContainer::ConstIterator PointDataIterator;

  /**  Type of the Transform Base class */
  typedef typename NumericTraits<CoordinateRepresentationType>::RealType TransformComputationType;
  typedef Transform<TransformComputationType,
                    itkGetStaticConstMacro(MovingMeshDimension),
                    itkGetStaticConstMacro(FixedMeshDimension)> TransformType;

  typedef typename TransformType::Pointer         TransformPointer;
  typedef typename TransformType::InputPointType  InputPointType;
  typedef typename TransformType::OutputPointType OutputPointType;
  typedef typename TransformType::ParametersType  TransformParametersType;
  typedef typename TransformType::JacobianType    TransformJacobianType;

  /**  Type of the Interpolator Base class */
  typedef InterpolateMeshFunction<MovingMeshType>   InterpolatorType;
  typedef typename InterpolatorType::Pointer        InterpolatorPointer;
  typedef typename InterpolatorType::RealType       RealDataType;
  typedef typename InterpolatorType::DerivativeType DerivativeDataType;

  /**  Type of the measure. */
  typedef Superclass::MeasureType MeasureType;

  /**  Type of the derivative. */
  typedef Superclass::DerivativeType DerivativeType;

  /**  Type of the parameters. */
  typedef Superclass::ParametersType ParametersType;

  /**  Type for the mask of the fixed image. Only pixels that are "inside"
       this mask will be considered for the computation of the metric */
  typedef SpatialObject<itkGetStaticConstMacro(FixedMeshDimension)> FixedMaskType;
  typedef typename  FixedMaskType::ConstPointer                     FixedMaskPointer;

  /**  Type for the mask of the moving image. Only pixels that are "inside"
       this mask will be considered for the computation of the metric */
  typedef SpatialObject<itkGetStaticConstMacro(MovingMeshDimension)> MovingMaskType;
  typedef typename  MovingMaskType::ConstPointer                     MovingMaskPointer;

  /** Connect the Fixed Pointset.  */
  itkSetConstObjectMacro( FixedMesh, FixedMeshType );

  /** Get the Fixed Pointset. */
  itkGetConstObjectMacro( FixedMesh, FixedMeshType );

  /** Connect the Moving Pointset.  */
  itkSetConstObjectMacro( MovingMesh, MovingMeshType );

  /** Get the Moving Pointset. */
  itkGetConstObjectMacro( MovingMesh, MovingMeshType );

  /** Connect the Transform. */
  itkSetObjectMacro( Transform, TransformType );

  /** Get a pointer to the Transform.  */
  itkGetObjectMacro( Transform, TransformType );

  /** Set the parameters defining the Transform. */
  void SetTransformParameters( const ParametersType & parameters ) const;

  /** Connect the Interpolator. */
  itkSetObjectMacro( Interpolator, InterpolatorType );

  /** Get a pointer to the Interpolator.  */
  itkGetConstObjectMacro( Interpolator, InterpolatorType );

  /** Return the number of parameters required by the Transform */
  unsigned int GetNumberOfParameters(void) const
  {
    return m_Transform->GetNumberOfParameters();
  }

  /** Connect the FixedMask .  */
  itkSetConstObjectMacro( FixedMask, FixedMaskType );

  /** Get the Fixed Pointset. */
  itkGetConstObjectMacro( FixedMask, FixedMaskType );

  /** Connect the Moving Pointset.  */
  itkSetConstObjectMacro( MovingMask, MovingMaskType );

  /** Get the Moving Pointset. */
  itkGetConstObjectMacro( MovingMask, MovingMaskType );

  /** Initialize the Metric by making sure that all the components
   *  are present and plugged together correctly     */
  virtual void Initialize(void) throw ( ExceptionObject );

protected:
  MeshToMeshMetric();
  virtual ~MeshToMeshMetric()
  {
  };
  void PrintSelf(std::ostream& os, Indent indent) const;

  FixedMeshConstPointer  m_FixedMesh;
  MovingMeshConstPointer m_MovingMesh;

  mutable TransformPointer m_Transform;
  InterpolatorPointer      m_Interpolator;

  mutable FixedMaskPointer  m_FixedMask;
  mutable MovingMaskPointer m_MovingMask;
private:
  MeshToMeshMetric(const Self &); // purposely not implemented
  void operator=(const Self &);   // purposely not implemented
};
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkMeshToMeshMetric.hxx"
#endif

#endif
