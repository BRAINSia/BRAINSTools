/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkMeshToMeshRegistrationMethod.h,v $
  Language:  C++
  Date:      $Date: 2006-06-08 20:30:47 $
  Version:   $Revision: 1.2 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkMeshToMeshRegistrationMethod_h
#define __itkMeshToMeshRegistrationMethod_h

#include "itkProcessObject.h"
#include "itkMeshToMeshMetric.h"
#include "itkSingleValuedNonLinearOptimizer.h"
#include "itkDataObjectDecorator.h"

namespace itk
{
/** \class MeshToMeshRegistrationMethod
 * \brief Base class for Mesh to Mesh Registration Methods
 *
 * This Class define the generic interface for a registration method.
 *
 * This class is templated over the type of the Mesh and the Mesh to be
 * registered. A generic Transform is used by this class. That allows to select
 * at run time the particular type of transformation that is to be applied for
 * registering the Meshs.
 *
 * This method use a generic Metric in order to compare the Mesh and the
 * Mesh.  The final goal of the registration method is to find the set of
 * parameters of the Transformation that optimizes the metric.
 *
 * The registration method also support a generic optimizer that can be
 * selected at run-time. The only restriction for the optimizer is that it
 * should be able to operate in single-valued cost functions given that the
 * metrics used to compare Mesh with Meshes provide a single value as
 * output.
 *
 * The terms : FixedMesh and MovingMesh are used in this class to indicate
 * that the Mesh is being mapped by the transform.
 *
 * This class uses the coordinate system of the Fixed Mesh as a reference
 * and searchs for a Transform that will map points from the space of the Fixed
 * Mesh to the space of the Moving Mesh.
 *
 * For doing so, a Metric will be continously applied to compare the Fixed
 * Mesh with the Transformed Moving Mesh. This process also requires to
 * interpolate values from the Moving Mesh.
 *
 * \ingroup RegistrationFilters
 */
template <typename TFixedMesh, typename TMovingMesh>
class ITK_EXPORT MeshToMeshRegistrationMethod : public ProcessObject
{
public:
  /** Standard class typedefs. */
  typedef MeshToMeshRegistrationMethod Self;
  typedef ProcessObject                Superclass;
  typedef SmartPointer<Self>           Pointer;
  typedef SmartPointer<const Self>     ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(MeshToMeshRegistrationMethod, ProcessObject);

  /**  Type of the Fixed Mesh. */
  typedef          TFixedMesh                  FixedMeshType;
  typedef typename FixedMeshType::ConstPointer FixedMeshConstPointer;

  /**  Type of the Moving Mesh. */
  typedef          TMovingMesh                  MovingMeshType;
  typedef typename MovingMeshType::ConstPointer MovingMeshConstPointer;

  /**  Type of the metric. */
  typedef MeshToMeshMetric<FixedMeshType,
                           MovingMeshType>   MetricType;
  typedef typename MetricType::Pointer MetricPointer;

  /**  Type of the Transform . */
  typedef  typename MetricType::TransformType TransformType;
  typedef  typename TransformType::Pointer    TransformPointer;

  /**  Type of the Interpolator. */
  typedef  typename MetricType::InterpolatorType InterpolatorType;
  typedef  typename InterpolatorType::Pointer    InterpolatorPointer;

  /** Type for the output: Using Decorator pattern for enabling
   *  the Transform to be passed in the data pipeline */
  typedef  DataObjectDecorator<TransformType>        TransformOutputType;
  typedef typename TransformOutputType::Pointer      TransformOutputPointer;
  typedef typename TransformOutputType::ConstPointer TransformOutputConstPointer;

  /**  Type of the optimizer. */
  typedef   SingleValuedNonLinearOptimizer OptimizerType;

  /** Type of the Transformation parameters This is the same type used to
   *  represent the search space of the optimization algorithm */
  typedef  typename MetricType::TransformParametersType ParametersType;

  /** Smart Pointer type to a DataObject. */
  typedef typename DataObject::Pointer DataObjectPointer;

  /** Method that initiates the registration. */
  void StartRegistration(void);

  /** Set/Get the Fixed Mesh. */
  itkSetConstObjectMacro( FixedMesh, FixedMeshType );
  itkGetConstObjectMacro( FixedMesh, FixedMeshType );

  /** Set/Get the Moving Mesh. */
  itkSetConstObjectMacro( MovingMesh, MovingMeshType );
  itkGetConstObjectMacro( MovingMesh, MovingMeshType );

  /** Set/Get the Optimizer. */
  itkSetObjectMacro( Optimizer,  OptimizerType );
  itkGetObjectMacro( Optimizer,  OptimizerType );

  /** Set/Get the Metric. */
  itkSetObjectMacro( Metric, MetricType );
  itkGetObjectMacro( Metric, MetricType );

  /** Set/Get the Transfrom. */
  itkSetObjectMacro( Transform, TransformType );
  itkGetObjectMacro( Transform, TransformType );

  /** Set/Get the Interpolator. */
  itkSetObjectMacro( Interpolator, InterpolatorType );
  itkGetObjectMacro( Interpolator, InterpolatorType );

  /** Set/Get the initial transformation parameters. */
  virtual void SetInitialTransformParameters( const ParametersType & param );

  itkGetConstReferenceMacro( InitialTransformParameters, ParametersType );

  /** Get the last transformation parameters visited by
   * the optimizer. */
  itkGetConstReferenceMacro( LastTransformParameters, ParametersType );

  /** Initialize by setting the interconnects between the components. */
  void Initialize() throw (ExceptionObject);

  /** Returns the transform resulting from the registration process  */
  const TransformOutputType * GetOutput() const;

  /** Make a DataObject of the correct type to be used as the specified
   * output. */
  virtual DataObjectPointer MakeOutput(size_t idx);

  /** Method to return the latest modified time of this object or
   * any of its cached ivars */
  unsigned long GetMTime() const;

protected:
  MeshToMeshRegistrationMethod();
  virtual ~MeshToMeshRegistrationMethod()
  {
  };
  void PrintSelf(std::ostream& os, Indent indent) const;

  /** Method invoked by the pipeline in order to trigger the computation of
   * the registration. */
  void  GenerateData();

private:
  MeshToMeshRegistrationMethod(const Self &); // purposely not implemented
  void operator=(const Self &);               // purposely not implemented

  MetricPointer          m_Metric;
  OptimizerType::Pointer m_Optimizer;

  MovingMeshConstPointer m_MovingMesh;
  FixedMeshConstPointer  m_FixedMesh;

  TransformPointer    m_Transform;
  InterpolatorPointer m_Interpolator;

  ParametersType m_InitialTransformParameters;
  ParametersType m_LastTransformParameters;
};
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkMeshToMeshRegistrationMethod.txx"
#endif

#endif
