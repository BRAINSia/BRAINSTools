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
#ifndef __itkPointLocator2_h
#define __itkPointLocator2_h

#include "itkObject.h"
#include "itkKdTree.h"
#include "itkKdTreeGenerator.h"
#include "itkPointSetToListSampleAdaptor.h"

namespace itk
{
/** \class PointLocator2
 * \brief Accelerate geometric searches for points.
 *
 * This class accelerates the search for the closest point to a user-provided
 * point, by using constructing a Kd-Tree structure for the PointSet.
 *
 */
template <class TPointSet>
class PointLocator2 : public Object
{
public:
  /** Standard class typedefs. */
  typedef PointLocator2            Self;
  typedef Object                   Superclass;
  typedef SmartPointer<Self>       Pointer;
  typedef SmartPointer<const Self> ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Standard part of every itk Object. */
  itkTypeMacro(PointLocator2, Object);

  itkStaticConstMacro(PointDimension, unsigned int, TPointSet::PointDimension);

  /** Typedefs related to the PointSet type */
  typedef TPointSet                                          PointSetType;
  typedef typename PointSetType::ConstPointer                PointSetConstPointer;
  typedef typename PointSetType::PointType                   PointType;
  typedef typename PointSetType::PointsContainer             PointsContainer;
  typedef typename PointSetType::PointsContainerPointer      PointsContainerPointer;
  typedef typename PointSetType::PointsContainerConstPointer PointsContainerConstPointer;

  typedef typename PointSetType::PointsContainerIterator      PointsContainerIteratorType;
  typedef typename PointSetType::PointsContainerConstIterator PointsContainerConstIteratorType;

  /** Define the specific PointSet Type */
  typedef itk::PointSet<typename TPointSet::PixelType, PointDimension> ListSamplePointSetType;

  /** Type of the PointSet to List Adaptor. */
  /** define the sample according to the point size */
  typedef itk::Statistics::PointSetToListSampleAdaptor<ListSamplePointSetType> SampleAdaptorType;

  typedef typename SampleAdaptorType::Pointer SampleAdaptorPointer;

  /** Types fo the KdTreeGenerator */
  typedef itk::Statistics::KdTreeGenerator<SampleAdaptorType> TreeGeneratorType;
  typedef typename TreeGeneratorType::Pointer                 TreeGeneratorPointer;
  typedef typename TreeGeneratorType::KdTreeType              TreeType;
  typedef typename TreeType::ConstPointer                     TreeConstPointer;
  typedef typename TreeType::InstanceIdentifierVectorType     InstanceIdentifierVectorType;

  /** Connect the PointSet as input */
  itkSetConstObjectMacro( PointSet, PointSetType );
  itkGetConstObjectMacro( PointSet, PointSetType );

  /** Pre-Compute the KdTree structure that will later facilitate the search of
   * points */
  void Initialize();

  /** Searches the k-nearest neighbors */
  void Search(const PointType & query, unsigned int numberOfNeighborsRequested,
              InstanceIdentifierVectorType& result) const;

  /** Searches the neighbors fallen into a hypersphere */
  void Search(const PointType & query, double radius, InstanceIdentifierVectorType& result) const;

protected:
  PointLocator2();
  ~PointLocator2();
  virtual void PrintSelf(std::ostream& os, Indent indent) const ITK_OVERRIDE;

private:
  ITK_DISALLOW_COPY_AND_ASSIGN(PointLocator2);

  PointSetConstPointer m_PointSet;
  SampleAdaptorPointer m_SampleAdaptor;
  TreeGeneratorPointer m_KdTreeGenerator;
  TreeConstPointer     m_Tree;
};
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkPointLocator2.hxx"
#endif

#endif
