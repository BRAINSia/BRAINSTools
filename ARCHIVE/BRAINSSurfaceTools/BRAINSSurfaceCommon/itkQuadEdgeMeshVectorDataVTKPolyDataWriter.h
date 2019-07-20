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
#ifndef __itkQuadEdgeMeshVectorDataVTKPolyDataWriter_h
#define __itkQuadEdgeMeshVectorDataVTKPolyDataWriter_h

#include "itkVTKPolyDataWriter.h"
#include <fstream>

namespace itk
{
/**
 * \class QuadEdgeMeshScalarDataVTKPolyData
 *
 * \brief This class saves a QuadMesh into a VTK-legacy file format,
 *        including its vector data associated with points.
 *
 * \ingroup Writers
 *
 */
template < typename TMesh >
class QuadEdgeMeshVectorDataVTKPolyDataWriter : public VTKPolyDataWriter< TMesh >
{
public:
  using Self = QuadEdgeMeshVectorDataVTKPolyDataWriter;
  using Pointer = SmartPointer< Self >;
  using ConstPointer = SmartPointer< const Self >;
  using Superclass = VTKPolyDataWriter< TMesh >;

  /** Run-time type information (and related methods).   */
  itkTypeMacro( QuadEdgeMeshVectorDataVTKPolyDataWriter, VTKPolyDataWriter );

  /** New macro for creation of through a Smart Pointer   */
  itkNewMacro( Self );

  using MeshType = TMesh;
  using MeshPointer = typename MeshType::Pointer;
  using CellType = typename MeshType::CellType;

  using PointsContainerPointer = typename MeshType::PointsContainerPointer;
  using PointsContainerIterator = typename MeshType::PointsContainerIterator;

  using PointDataContainerPointer = typename MeshType::PointDataContainerPointer;
  using PointDataContainerConstPointer = typename MeshType::PointDataContainerConstPointer;
  using PointDataContainerIterator = typename MeshType::PointDataContainerIterator;

  using CellsContainer = typename MeshType::CellsContainer;
  using CellsContainerPointer = typename CellsContainer::Pointer;
  using CellsContainerConstPointer = typename CellsContainer::ConstPointer;
  using CellsContainerIterator = typename CellsContainer::Iterator;
  using CellsContainerConstIterator = typename CellsContainer::ConstIterator;

  using CellDataContainer = typename MeshType::CellDataContainer;
  using CellDataContainerPointer = typename CellDataContainer::Pointer;
  using CellDataContainerConstPointer = typename CellDataContainer::ConstPointer;
  using CellDataContainerIterator = typename CellDataContainer::Iterator;
  using CellDataContainerConstIterator = typename CellDataContainer::ConstIterator;

  /** Set/Get the name of the CellDataName where data are written. */
  itkSetStringMacro( CellDataName );
  itkGetStringMacro( CellDataName );

  /** Set/Get the name of the PointDataName where data are written. */
  itkSetStringMacro( PointDataName );
  itkGetStringMacro( PointDataName );

protected:
  QuadEdgeMeshVectorDataVTKPolyDataWriter();
  ~QuadEdgeMeshVectorDataVTKPolyDataWriter();

  std::string m_CellDataName;
  std::string m_PointDataName;

  void
  GenerateData() override;

  void
  WriteCellData();

  void
  WritePointData();

private:
  QuadEdgeMeshVectorDataVTKPolyDataWriter( const Self & );
  void
  operator=( const Self & );
};
} // namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#  include "itkQuadEdgeMeshVectorDataVTKPolyDataWriter.hxx"
#endif

#endif
