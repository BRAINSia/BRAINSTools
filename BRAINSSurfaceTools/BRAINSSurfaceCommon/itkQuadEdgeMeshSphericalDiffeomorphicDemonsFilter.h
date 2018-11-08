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
#ifndef __itkQuadEdgeMeshSphericalDiffeomorphicDemonsFilter_h
#define __itkQuadEdgeMeshSphericalDiffeomorphicDemonsFilter_h

#include "itkQuadEdgeMeshToQuadEdgeMeshFilter.h"
#include "itkNodeScalarGradientCalculator.h"
#include "itkNodeVectorJacobianCalculator.h"
#include "itkLinearInterpolateDeformationFieldMeshFunction.h"
#include "itkTriangleListBasisSystemCalculator.h"
#include "itkTriangleBasisSystem.h"
#include "itkVectorContainer.h"
#include "itkVector.h"
#include "itkTimeProbesCollectorBase.h"

namespace itk
{
template <typename TFixedMesh, typename TMovingMesh, typename TOutputMesh>
class QuadEdgeMeshSphericalDiffeomorphicDemonsFilter :
  public QuadEdgeMeshToQuadEdgeMeshFilter<TFixedMesh, TOutputMesh>
{
public:
  typedef QuadEdgeMeshSphericalDiffeomorphicDemonsFilter            Self;
  typedef SmartPointer<Self>                                        Pointer;
  typedef SmartPointer<const Self>                                  ConstPointer;
  typedef QuadEdgeMeshToQuadEdgeMeshFilter<TFixedMesh, TOutputMesh> Superclass;

  /** Method that instantiates a new object */
  itkNewMacro( Self );

  /** Method that provides the name of the class as a string as well as the
   * name of the parent class. */
  itkTypeMacro( QuadEdgeMeshSphericalDiffeomorphicDemonsFilter, QuadEdgeMeshToQuadEdgeMeshFilter );

  /** Input types. */
  typedef TFixedMesh                                       FixedMeshType;
  typedef TMovingMesh                                      MovingMeshType;
  typedef typename  FixedMeshType::Pointer                 FixedMeshPointer;
  typedef typename  FixedMeshType::ConstPointer            FixedMeshConstPointer;
  typedef typename  FixedMeshType::PointType               PointType;
  typedef typename  FixedMeshType::PixelType               FixedPixelType;
  typedef typename  FixedMeshType::PointsContainer         FixedPointsContainer;
  typedef typename  FixedPointsContainer::Iterator         FixedPointsIterator;
  typedef typename  FixedPointsContainer::ConstIterator    FixedPointsConstIterator;
  typedef typename  FixedMeshType::PointDataContainer      FixedPointDataContainer;
  typedef typename  FixedPointDataContainer::ConstIterator FixedPointDataConstIterator;
  typedef typename  MovingMeshType::ConstPointer           MovingMeshConstPointer;

  static constexpr unsigned int PointDimension = FixedMeshType::PointDimension;

  /** Output types. */
  typedef TOutputMesh                                    OutputMeshType;
  typedef typename  OutputMeshType::Pointer              OutputMeshPointer;
  typedef typename  Superclass::OutputPointDataContainer OutputPointDataContainer;
  typedef typename OutputPointDataContainer::Pointer     OutputPointDataContainerPointer;
  typedef typename OutputPointDataContainer::Iterator    OutputPointDataContainerIterator;

  /** Declaration of internal types, some of which are exposed for monitoring purposes */
  typedef typename PointType::VectorType                    VectorType;
  typedef TriangleBasisSystem<VectorType, 2>                BasisSystemType;
  typedef typename FixedMeshType::PointIdentifier           PointIdentifier;
  typedef VectorContainer<PointIdentifier, BasisSystemType> BasisSystemContainerType;
  typedef typename BasisSystemContainerType::Pointer        BasisSystemContainerPointer;
  typedef typename BasisSystemContainerType::Iterator       BasisSystemContainerIterator;
  typedef typename FixedMeshType::Traits                    FixedMeshTraits;
  typedef PointSet<FixedPixelType, PointDimension, FixedMeshTraits>
    DestinationPointSetType;
  typedef typename DestinationPointSetType::PointsContainer     DestinationPointContainerType;
  typedef typename DestinationPointContainerType::Pointer       DestinationPointContainerPointer;
  typedef typename DestinationPointContainerType::Iterator      DestinationPointIterator;
  typedef typename DestinationPointContainerType::ConstIterator DestinationPointConstIterator;
  typedef VectorContainer<PointIdentifier, double>              NodeSigmaContainerType;
  typedef typename NodeSigmaContainerType::Pointer              NodeSigmaContainerPointer;
  typedef typename NodeSigmaContainerType::Iterator             NodeSigmaContainerIterator;
  typedef typename NodeSigmaContainerType::ConstPointer         NodeSigmaContainerConstPointer;
  typedef typename NodeSigmaContainerType::ConstIterator        NodeSigmaContainerConstIterator;

  typedef Vector<double, 3>                                    VelocityVectorType;
  typedef VectorContainer<PointIdentifier, VelocityVectorType> VelocityVectorContainer;
  typedef typename VelocityVectorContainer::Pointer            VelocityVectorPointer;
  typedef typename VelocityVectorContainer::Iterator           VelocityVectorIterator;
  typedef typename VelocityVectorContainer::ConstPointer       VelocityVectorConstPointer;
  typedef typename VelocityVectorContainer::ConstIterator      VelocityVectorConstIterator;

  typedef Vector<double, 3>                                   TangentVectorType;
  typedef VectorContainer<PointIdentifier, TangentVectorType> TangentVectorContainer;
  typedef typename TangentVectorContainer::Pointer            TangentVectorPointer;
  typedef typename TangentVectorContainer::Iterator           TangentVectorIterator;
  typedef typename TangentVectorContainer::ConstPointer       TangentVectorConstPointer;
  typedef typename TangentVectorContainer::ConstIterator      TangentVectorConstIterator;

  typedef VectorContainer<PointIdentifier, double>            ShortestLengthContainerType;
  typedef typename ShortestLengthContainerType::Pointer       ShortestLengthContainerPointer;
  typedef typename ShortestLengthContainerType::Iterator      ShortestLengthContainerIterator;
  typedef typename ShortestLengthContainerType::ConstIterator ShortestLengthContainerConstIterator;

  /** Set/Get the Fixed mesh. */
  void SetFixedMesh( const FixedMeshType * fixedMesh );

  itkGetConstObjectMacro( FixedMesh, FixedMeshType );

  /** Set/Get the Moving mesh. */
  void SetMovingMesh( const MovingMeshType * movingMesh );

  itkGetConstObjectMacro( MovingMesh, MovingMeshType );

  /** Returns the array of local coordinates systems at every node of the fixed
   * mesh. This array is only valid after a call to Update() has completed
   * successfully. */
  itkGetConstObjectMacro( BasisSystemAtNode, BasisSystemContainerType );

  /** Returns the array of destination points used for initializing the
   * deformation field to all nodes of the Fixed Mesh. The points are stored
   * in the point container of an itk::PointSet. */
  const DestinationPointSetType * GetInitialDestinationPoints() const;

  /** Set the initial list of destination points. This makes possible to use the
   * outcome of a deformable registration as the input for another one. It is
   * also possible to initialize such array of destination points by using a
   * low degrees of freedom transform, such as a rigid transform, for example.
   * The destination points are stored as the points of an itk::PointSet class.
   * This is used as a data pipeline input, meaning that when calling Update()
   * in this filter, the pipeline is going to verify if the content of the
   * PointSet is computed by a preceding filter, and in such case, whether
   * the preceding filter must be re-executed or not. */
  void SetInitialDestinationPoints( const DestinationPointSetType * );

  /** Returns the array of destination points resulting from applying the
   * deformation field to all nodes of the Fixed Mesh. The points are stored
   * in the point container of an itk::PointSet. The content of this array
   * is only valid after the first iteration of the filter execution has been
   * completed. It can be used for tracking the progress of the filter. The
   * returned point is non-const just to allow users to call
   * DisconnectPipeline(). */
  DestinationPointSetType * GetFinalDestinationPoints() const;

  /** Set/Get the maximum number of iterations that the filter will be
   * allowed to run.  The default is set to 50. */
  itkSetMacro( MaximumNumberOfIterations, unsigned int );
  itkGetMacro( MaximumNumberOfIterations, unsigned int );

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

  /** Set/Get the value of the regularization constant used in the computation
   * of the deformation field update. */
  itkSetMacro( Epsilon, double );
  itkGetConstMacro( Epsilon, double );

  /** Set/Get the value of the weight used in the contribution of the Jacobian
   * to the Levenberg Marquardt term during the computation of the velocity
   * field. */
  itkSetMacro( SigmaX, double );
  itkGetConstMacro( SigmaX, double );

  /** Variable that defines whether the filter will self-adjust the values of
   * SigmaX and Epsilon in order to get closer to the ratio of
   * largestVelocityMagnitude being similar to the value of the shortest edge
   * length. Default: false*/
  itkSetMacro( SelfRegulatedMode, bool );
  itkGetConstMacro( SelfRegulatedMode, bool );
  itkBooleanMacro( SelfRegulatedMode );

  /** Variable that defines whether the filter will stop the iterations
  by checking the change of metric between the previous and current
  iterations.*/
  itkSetMacro( SelfStopMode, bool );
  itkGetConstMacro( SelfStopMode, bool );
  itkBooleanMacro( SelfStopMode );

  /** Set/Get the container of sigma values to be associated with each node of
   * the fixed mesh. This sigma value represents the expected variability of
   * scalar values at this node of the mesh. */
  itkSetConstObjectMacro( FixedNodesSigmas, NodeSigmaContainerType );
  itkGetConstObjectMacro( FixedNodesSigmas, NodeSigmaContainerType );

  /** The smoothing filter will run iteratively until reaching this maximum
   * number of iterations. Emprical observartions indicate that ten iterations
   * are enough for typical deformation fields, but of course this would depend
   * on the process that you used for generating your deformation field.
   */
  itkSetMacro( MaximumNumberOfSmoothingIterations, unsigned int );
  itkGetConstMacro( MaximumNumberOfSmoothingIterations, unsigned int );

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
  itkGetConstMacro( Lambda, double );

  /** Returns a new mesh that is identical to the input Fixed Mesh except for
   * the fact that its points have been moved to the location of the destination
   * points. This is equivalent to deforming the input fixed mesh by applying to
   * it the deformation field that is computed by this filter.  This output is
   * only valid after you have successfully called the Update() method. */
  const FixedMeshType * GetDeformedFixedMesh() const;

  /** Square root of sum of squared differences of intensities between nodes in
   * the mesh. This value is not used directly in the algorithm. We simply
   * compute it as a way of providing feedback on the progress of the
   * registration. */
  itkGetConstMacro( MetricValue, double );

  /** m_MetricChange =
  100 * (m_MetricValue(n) - m_MetricValue(n-1))/m_MetricValue(n-1)
  n: current iteration; n-1: previous iteration*/
  itkGetConstMacro( MetricChange, double );

  /** The significant difference between current iteration and the previous iteration.
  If m_MetricChange is bigger than it, iterations keep going until it hits the
  m_MaximumNumberOfIterations. If m_MetricChange is smaller than it, registration
  stops.*/
  itkSetMacro( MetricSignificance, double );
  itkGetConstMacro( MetricSignificance, double );

  using Superclass::MakeOutput;
  /**  Create the Output of the proper type for that output number */
  virtual DataObject::Pointer MakeOutput(size_t idx) override;

  /** Print out in the argument ostream the results of the chronometer
   * measurements. This is intended to be used for profiling the deformation
   * filter. */
  void ChronometerReport( std::ofstream & os ) const;

protected:
  QuadEdgeMeshSphericalDiffeomorphicDemonsFilter();
  ~QuadEdgeMeshSphericalDiffeomorphicDemonsFilter();
  void PrintSelf(std::ostream& os, Indent indent) const override;

  virtual void GenerateData() override;

private:
  QuadEdgeMeshSphericalDiffeomorphicDemonsFilter( const Self & );
  void operator =( const Self & );

  void AllocateInternalArrays();

  void InitializeFixedNodesSigmas();

  void ComputeBasisSystemAtEveryNode();

  void ComputeInitialArrayOfDestinationPoints();

  void InitializeInterpolators();

  void InitializeGradientCalculators();

  void CopyInitialDestinationPoints();

  void CopySourcePoinstAsDestinationPoints();

  void RunIterations();

  void ComputeMappedMovingValueAtEveryNode();

  void ComputeGradientsOfMappedMovingValueAtEveryNode();

  void ComputeVelocityField();

  void ComputeSelfRegulatedVelocityField();

  void SmoothDeformationField();

  void ConvertDeformationFieldToTangentVectorField();

  void SmoothTangentVectorField();

  void ConvertTangentVectorFieldToDeformationField();

  void AssignResampledMovingValuesToOutputMesh();

  void ComposeFixedMeshOutputDisplacedToMovingMesh();

  void ComposeDestinationPointsOutputPointSet();

  void CopyDestinationPointsToDeformedFixedMesh();

  void ComputeScalingAndSquaringNumberOfIterations();

  void ComputeShortestEdgeLength();

  double ComputeLargestVelocityMagnitude() const;

  void ComputeLargestVelocityMagnitudeToShortestEdgeLengthRatio();

  void ComputeDeformationByScalingAndSquaring();

  void ComposeDeformationUpdateWithPreviousDeformation();

  void ComputeSelfRegulatedSigmaXandEpsilon();

  void SwapOldAndNewDestinationPointContainers();

  void SwapOldAndNewDisplacementFieldContainers();

  void SwapOldAndNewTangetFieldContainers();

  void ParalelTransport(const PointType sourcePoint, const PointType destinationPoint,
                        const TangentVectorType & inputVector, TangentVectorType & transportedVector ) const;

  void PrintOutDeformationVectors( std::ostream & os = std::cout );

  virtual PointType InterpolateDestinationFieldAtPoint(const DestinationPointContainerType * destinationField,
                                                       const PointType & point );

  virtual void ProjectPointToSphereSurface( PointType & point ) const;

  MovingMeshConstPointer m_MovingMesh;
  FixedMeshConstPointer  m_FixedMesh;
  FixedMeshPointer       m_FixedMeshAtInitialDestinationPoints;

  /** This is the Array of "Qn" matrices
   *  presented in equation 3.14 in the paper. */
  BasisSystemContainerPointer m_BasisSystemAtNode;

  /** Array containing the destination coordinates of every node in the Fixed
   * Mesh.  This array represents both the deformation field c(xn) and its
   * smoothed version, the field s(xn) as defined in.  */
  DestinationPointContainerPointer m_DestinationPoints;
  DestinationPointContainerPointer m_DestinationPointsSwap;
  bool                             m_UserProvidedInitialDestinationPoints;

  /** Auxiliary array for computing the Exponential of the velocity field
   * via the Scaling and Squaring method.  */
  DestinationPointContainerPointer m_DisplacementField;
  DestinationPointContainerPointer m_DisplacementFieldSwap;

  /** Maximum number of iterations that the filter will be allowed to run. */
  unsigned int m_MaximumNumberOfIterations;

  /** Maximum number of iterations that will be used for smoothing the tangent field. */
  unsigned int m_MaximumNumberOfSmoothingIterations;

  /** Coefficient that controls the degree of smoothing applied to the tangent field. */
  double m_Lambda;

  typedef TriangleListBasisSystemCalculator<
      FixedMeshType, BasisSystemType> TriangleListBasisSystemCalculatorType;

  /** Helper class that will compute basis systems at every triangle of the Fixed Mesh. */
  typename TriangleListBasisSystemCalculatorType::Pointer m_TriangleListBasisSystemCalculator;

  /** Types definitions for the container of values resampled from the Moving
   * mesh into the coordinates of the Fixed mesh nodes. */
  typedef typename MovingMeshType::PixelType                    MovingPixelType;
  typedef typename NumericTraits<MovingPixelType>::RealType     MovingPixelRealType;
  typedef VectorContainer<PointIdentifier, MovingPixelRealType> ResampledMovingValuesContainerType;
  typedef typename ResampledMovingValuesContainerType::Iterator ResampledMovingValuesContainerIterator;

  typedef typename NumericTraits<FixedPixelType>::RealType FixedPixelRealType;

  /** Container that stores values resampled from the Moving mesh field at the
   * coordinates resulting from mapping the fixed mesh nodes through the current
   * deformation field. */
  typename ResampledMovingValuesContainerType::Pointer          m_ResampledMovingValuesContainer;

  /** Interpolator type for bringing scalar values from the Moving Mesh into the Fixed Mesh. */
  typedef LinearInterpolateMeshFunction<MovingMeshType> ScalarInterpolatorType;

  /** Interpolator object that will bring scalar values from the Moving Mesh into the Fixed Mesh. */
  typename ScalarInterpolatorType::Pointer                      m_ScalarInterpolator;

  /** Interpolator for the deformation field values on the grid of the Fixed mesh. */
  typedef LinearInterpolateDeformationFieldMeshFunction<FixedMeshType> DeformationInterpolatorType;

  /** Interpolator object that will compute deformation destination points on the fixed mesh grid. */
  typename DeformationInterpolatorType::Pointer                 m_DeformationInterpolator;

  /** Helper class that will compute the gradient of resampled Moving mesh
   * values at every node of the Fixed mesh with respect to the coordinate system
   * of that node in the fixed mesh. */
  typedef NodeScalarGradientCalculator<
      FixedMeshType, ResampledMovingValuesContainerType>         NodeScalarGradientCalculatorType;
  typename NodeScalarGradientCalculatorType::Pointer            m_NodeScalarGradientCalculator;

  /** Helper class that will compute the Jacobian of destination points
   * at every node of the Fixed mesh with respect to the coordinate system
   * of that node in the fixed mesh. */
  typedef NodeVectorJacobianCalculator<
      FixedMeshType, DestinationPointContainerType>              NodeVectorJacobianCalculatorType;
  typename NodeVectorJacobianCalculatorType::Pointer            m_NodeVectorJacobianCalculator;

  /** Center of spherical mesh. We assume that both the Fixed and
   * Moving meshes have spherical geometry and that they share the same
   * center and radius. */
  PointType m_SphereCenter;

  /** Radius of spherical mesh. We assume that both the Fixed and
   * Moving meshes have spherical geometry and that they share the same
   * center and radius. */
  double m_SphereRadius;

  /** Regularization constant used during the update of the deformation field. */
  double m_Epsilon;          // 1/sigmaX^2 == epsilon === > control step length of optimization.

  /** This term controls the contribution of the Jacobian in the Levenberg
   * Marquardt term that computes the velocity field. Large values of SigmaX
   * will result in large deformations. */
  double m_SigmaX;

  /** Shortest edge length in the input fixed mesh. This is used for computing
   * the number of iterations required in the squaring and scaling method. */
  double m_ShortestEdgeLength;

  /** Number of iterations to be applied in the Scaling and Squaring method. */
  unsigned int m_ScalingAndSquaringNumberOfIterations;

  /** Container of Sigma values associated to each node of the Fixed mesh. */
  NodeSigmaContainerConstPointer m_FixedNodesSigmas;

  /** Container of velocity field vectors. */
  VelocityVectorPointer m_VelocityField;

  /** Container of tangent vectors used to smooth the deformation field. */
  TangentVectorPointer m_TangentVectorField;
  TangentVectorPointer m_TangentVectorFieldSwap;

  TimeProbesCollectorBase m_Chronometer;

  /** Square root of sum of squared differences of intensities between nodes in
   * the mesh. This value is not used directly in the algorithm. We simply
   * compute it as a way of providing feedback on the progress of the
   * registration. */
  double m_MetricValue;

  /** m_MetricChange =
  100 * (m_MetricValue(n) - m_MetricValue(n-1))/m_MetricValue(n-1)
  n: current iteration; n-1: previous iteration*/
  double m_MetricChange;

  /** The significant difference between current iteration and the previous iteration.
  If m_MetricChange is bigger than it, iterations keep going until it hits the
  m_MaximumNumberOfIterations. If m_MetricChange is smaller than it, registration
  stops.*/
  double m_MetricSignificance;

  /** Variable that defines whether the filter will self-adjust the values of
   * SigmaX and Epsilon in order to get closer to the ratio of
   * largestVelocityMagnitude being similar to the value of the shortest edge
   * length. */
  bool m_SelfRegulatedMode;

  /** Variable that defines whether the filter will stop the iterations
  by checking the change of metric between the previous and current
  iterations.*/
  bool m_SelfStopMode;

  /** largest ratio of velocity versus shortest edge length of the
   * corresponding node. */
  double m_LargestVelocityToEdgeLengthRatio;

  /** Container of lengths corresponding to the shortest edge of every node. */
  ShortestLengthContainerPointer m_ShortestEdgeLengthPerPoint;
};
}

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkQuadEdgeMeshSphericalDiffeomorphicDemonsFilter.hxx"
#endif

#endif
