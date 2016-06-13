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
#ifndef __itkNodeScalarGradientCalculator_h
#define __itkNodeScalarGradientCalculator_h

#include "itkFunctionBase.h"
#include "itkCovariantVector.h"
#include "itkTriangleBasisSystem.h"
#include "itkTriangleListBasisSystemCalculator.h"
#include "itkLinearInterpolateMeshFunction.h"
#include "itkTriangleHelper.h"

namespace itk
{
/** \class NodeScalarGradientCalculator
 * Estimates Gradient of a scalar function at at the node of a mesh.
 *
 * \brief
 *
 * \sa NodeScalarGradientCalculator
 * \ingroup MeshFunctions
 *
 * */
template <class TInputMesh, class TPointDataContainer = double>
class NodeScalarGradientCalculator :
  public         FunctionBase<typename TInputMesh::PointIdentifier,
                              CovariantVector<
                                typename NumericTraits<typename TPointDataContainer::Element>::RealType,
                                TInputMesh::PointDimension> >
{
public:
  /** Standard class typedefs. */
  typedef NodeScalarGradientCalculator Self;
  typedef FunctionBase<typename TInputMesh::PointIdentifier,
                       CovariantVector<
                         typename NumericTraits<typename TPointDataContainer::Element>::RealType,
                         TInputMesh::PointDimension> >  Superclass;
  typedef SmartPointer<Self>       Pointer;
  typedef SmartPointer<const Self> ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(NodeScalarGradientCalculator, FunctionBase);

  /** Dimension underlying input mesh. */
  itkStaticConstMacro(MeshDimension, unsigned int,
                      TInputMesh::PointDimension );

  /** Point typedef support. */
  typedef TInputMesh                                 InputMeshType;
  typedef typename InputMeshType::PointType          PointType;
  typedef typename PointType::VectorType             VectorType;
  typedef typename InputMeshType::PointDataContainer PointDataContainer;
  typedef typename PointDataContainer::Iterator      PointDataIterator;
  typedef typename InputMeshType::PixelType          PixelType;
  typedef typename InputMeshType::PointsContainer    PointsContainer;
  typedef typename PointsContainer::ConstIterator    PointIterator;

  typedef TriangleBasisSystem<VectorType, 2>                       TriangleBasisSystemType;
  typedef typename InputMeshType::CellIdentifier                   CellIdentifier;
  typedef VectorContainer<CellIdentifier, TriangleBasisSystemType> BasisSystemListType;
  typedef typename BasisSystemListType::ConstIterator              BasisSystemListIterator;

  typedef typename InputMeshType::CellType       CellType;
  typedef typename InputMeshType::CellTraits     CellTraits;
  typedef typename InputMeshType::CellsContainer CellsContainer;
  typedef typename CellsContainer::Iterator      CellsContainerIterator;
  typedef typename CellsContainer::ConstIterator CellsContainerConstIterator;
  typedef typename CellTraits::PointIdIterator   PointIdIterator;

  typedef TriangleHelper<PointType>                 TriangleType;
  typedef typename TriangleType::CoordRepType       AreaType;
  typedef VectorContainer<CellIdentifier, AreaType> AreaListType;
  typedef typename AreaListType::Iterator           AreaListIterator;
  typedef typename AreaListType::ConstIterator      AreaListConstIterator;

  typedef TriangleListBasisSystemCalculator<InputMeshType, TriangleBasisSystemType>
    TriangleListBasisSystemCalculatorType;

  typedef LinearInterpolateMeshFunction<InputMeshType> InterpolatorType;

  typedef PointLocator2<TInputMesh>                  PointLocatorType;
  typedef typename PointLocatorType::Pointer         PointLocatorPointer;
  typedef typename InterpolatorType::PointIdentifier PointIdentifier;

  typedef typename InterpolatorType::RealType              RealType;
  typedef typename InterpolatorType::DerivativeType        DerivativeType;
  typedef VectorContainer<PointIdentifier, DerivativeType> DerivativeListType;
  typedef typename DerivativeListType::Iterator            DerivativeListIterator;
  typedef typename DerivativeListType::ConstIterator       DerivativeListConstIterator;

  typedef typename PointType::CoordRepType               CoordRepType;
  typedef VectorContainer<PointIdentifier, CoordRepType> CoordRepListType;

  /** Set/Get the input mesh. */
  itkSetConstObjectMacro( InputMesh, InputMeshType );
  itkGetConstObjectMacro( InputMesh, InputMeshType );

  /** Set/Get the input mesh. */
  itkSetConstObjectMacro( DataContainer, TPointDataContainer );
  itkGetConstObjectMacro( DataContainer, TPointDataContainer );

  /** Definition of input type and output type. The input type is actually a
   * point identifier, while the output type is a gradient of the scalar values. */
  typedef typename Superclass::OutputType OutputType;
  typedef typename Superclass::InputType  InputType;

  /** Initialize internal variables. This method MUST be called first, before
   * calling the Compute() method. In a normal usage, the Initialize() method
   * will be called once ant the beginning, and then it will be followed by
   * multiple calls to the Compute() method.
   *
   * \sa Compute
   */
  virtual void Initialize( void );

  /** Compute the values at all the nodes. The Initialize() method MUST be
   * called first. In a normal usage, the Initialize() method will be called
   * once ant the beginning, and then it will be followed by multiple calls to
   * the Compute() method.
   *
   * \sa Initialize
   */
  virtual void Compute();

  /** Set/Get list of basis systems at every cell. */
  itkSetConstObjectMacro( BasisSystemList, BasisSystemListType );
  itkGetConstObjectMacro( BasisSystemList, BasisSystemListType );

  /** Evaluate at the specified input position */
  virtual OutputType Evaluate( const InputType& input) const ITK_OVERRIDE;

  /** Set Sphere Center.  The implementation of this class assumes that the
   * Mesh surface has a spherical geometry (not only spherical topology). With
   * this method you can specify the coordinates of the center of the sphere
   * represented by the Mesh. This will be used in the computation of parallel
   * transport for vector values associated with nodes.
   */
  itkSetMacro( SphereCenter, PointType );
  itkGetConstMacro( SphereCenter, PointType );

  /** Set Sphere Radius.  The implementation of this class assumes that the
   * Mesh surface has a spherical geometry (not only spherical topology). With
   * this method you can specify the radius of the sphere. This will be used in
   * the computation of parallel transport for vector values associated
   * with nodes.
   */
  itkSetMacro( SphereRadius, double );
  itkGetConstMacro( SphereRadius, double );
protected:
  NodeScalarGradientCalculator();
  ~NodeScalarGradientCalculator();

  void PrintSelf(std::ostream& os, Indent indent) const ITK_OVERRIDE;

private:
  ITK_DISALLOW_COPY_AND_ASSIGN(NodeScalarGradientCalculator);

  typename InputMeshType::ConstPointer                 m_InputMesh;
  typename TPointDataContainer::ConstPointer           m_DataContainer;
  typename BasisSystemListType::ConstPointer           m_BasisSystemList;
  typename AreaListType::Pointer                       m_AreaList;
  typename CoordRepListType::Pointer                   m_PointAreaAccumulatorList;
  typename DerivativeListType::Pointer                 m_PointDerivativeAccumulatorList;

  /** Check that all necessary inputs are connected. */
  virtual void VerifyInputs( void ) const;

  /** Allocate memory for all the internal containers. This method is called
   * only once during the initialization of the calculator. It doesn't need to be
   * called again unless the number of cells or number of points in the input
   * mesh change. */
  void AllocateInternalContainers();

  /** Parallel-transport for gradient vectors. This is equivalent to sliding them along
   * the great circle that connects sourcePoint with destinationPoint.  */
  void ParalelTransport( const PointType sourcePoint, const PointType destinationPoint,
                         const DerivativeType & inputVector, DerivativeType & transportedVector ) const;

  /** Divide the cumulated derivatives by the cumulated areas */
  void NormalizeDerivativesByTotalArea();

  /** Compute the area in all cells and store it in a container */
  void ComputeAreaForAllCells();

  /** Fill the values of several containers with null values. This is done as
   * initialization before we start accumulating values in them. */
  void SetContainersToNullValues();

  PointType m_SphereCenter;
  double    m_SphereRadius;
};
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkNodeScalarGradientCalculator.hxx"
#endif

#endif
