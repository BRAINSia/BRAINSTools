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
#ifndef __itkPointLocator2_hxx
#define __itkPointLocator2_hxx
#include "itkPointLocator2.h"

namespace itk
{
template <class TMesh>
PointLocator2<TMesh>
::PointLocator2()
{
  this->m_SampleAdaptor = SampleAdaptorType::New();
  this->m_KdTreeGenerator = TreeGeneratorType::New();
}

template <class TMesh>
PointLocator2<TMesh>
::~PointLocator2()
{
}

template <class TMesh>
void
PointLocator2<TMesh>
::Initialize()
{
  this->m_SampleAdaptor = SampleAdaptorType::New();
  this->m_KdTreeGenerator = TreeGeneratorType::New();

  // Lack of const-correctness in the PointSetAdaptor should be fixed.
  typename ListSamplePointSetType::Pointer pointSet = ListSamplePointSetType::New();

  PointsContainerConstPointer pointsContainer = this->m_PointSet->GetPoints();

  PointsContainerConstIteratorType pIr = pointsContainer->Begin();
  PointsContainerConstIteratorType pIrEnd = pointsContainer->End();

  typename ListSamplePointSetType::PointType point;
  int i = 0;
  while( pIr != pIrEnd )
    {
    point.CastFrom( pIr.Value() );
    pointSet->SetPoint(i, point);
    i++; pIr++;
    }

  this->m_SampleAdaptor->SetPointSet(pointSet.GetPointer() );
  // this->m_SampleAdaptor->SetPointSet(
  //  const_cast< PointSetType * >( this->m_PointSet.GetPointer() ) );

  // this->m_SampleAdaptor->SetMeasurementVectorSize( PointDimension );

  this->m_KdTreeGenerator->SetSample( this->m_SampleAdaptor );
  this->m_KdTreeGenerator->SetBucketSize( 16 );

  this->m_KdTreeGenerator->Update();

  this->m_Tree = this->m_KdTreeGenerator->GetOutput();
}

template <class TMesh>
void
PointLocator2<TMesh>
::Search(const PointType & query,
         unsigned int numberOfNeighborsRequested,
         InstanceIdentifierVectorType& result) const
{
  this->m_Tree->Search( query, numberOfNeighborsRequested, result );
}

template <class TMesh>
void
PointLocator2<TMesh>
::Search(const PointType & query,
         double radius,
         InstanceIdentifierVectorType& result) const
{
  this->m_Tree->Search( query, radius, result );
}

/**
 * Print out internals
 */
template <class TMesh>
void
PointLocator2<TMesh>
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);
}
} // end namespace itk

#endif
