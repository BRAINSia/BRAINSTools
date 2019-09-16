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
#ifndef __itkQuadEdgeMeshVTKPolyDataReader_h
#define __itkQuadEdgeMeshVTKPolyDataReader_h

#include "itkMesh.h"
#include "itkMeshSource.h"
#include "itkTriangleCell.h"
#include "itkMapContainer.h"

namespace itk
{
/** \class QuadEdgeMeshVTKPolyDataReader
 * \brief
 * Reads a vtkPolyData file and create an itkMesh.
 *
 * Caveat: itkQuadEdgeMeshVTKPolyDataReader can only read triangle meshes.
 *         Use vtkTriangleFilter to convert your mesh to a triangle mesh.
 */
template <typename TOutputMesh>
class QuadEdgeMeshVTKPolyDataReader : public MeshSource<TOutputMesh>
{
public:
  ITK_DISALLOW_COPY_AND_ASSIGN(QuadEdgeMeshVTKPolyDataReader);

  /** Standard "Self" type alias. */
  using Self = QuadEdgeMeshVTKPolyDataReader;
  using Superclass = MeshSource<TOutputMesh>;
  using Pointer = SmartPointer<Self>;
  using ConstPointer = SmartPointer<const Self>;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(QuadEdgeMeshVTKPolyDataReader, MeshSource);

  /** Hold on to the type information specified by the template parameters. */
  using OutputMeshType = TOutputMesh;
  using MeshTraits = typename OutputMeshType::MeshTraits;
  using PointType = typename OutputMeshType::PointType;
  using PixelType = typename MeshTraits::PixelType;

  /** Some convenient type alias. */
  using OutputMeshPointer = typename OutputMeshType::Pointer;
  using CellTraits = typename OutputMeshType::CellTraits;
  using CellIdentifier = typename OutputMeshType::CellIdentifier;
  using CellType = typename OutputMeshType::CellType;
  using CellAutoPointer = typename OutputMeshType::CellAutoPointer;
  using PointIdentifier = typename OutputMeshType::PointIdentifier;
  using PointIdIterator = typename CellTraits::PointIdIterator;

  typedef typename OutputMeshType::PointsContainerPointer PointsContainerPointer;

  typedef typename OutputMeshType::PointsContainer PointsContainer;

  /** Define the triangular cell types which form the surface  */
  using TriangleCellType = TriangleCell<CellType>;

  typedef typename TriangleCellType::SelfAutoPointer TriangleCellAutoPointer;

  using IndexPairType = std::pair<unsigned long, unsigned long>;
  using PointMapType = MapContainer<IndexPairType, unsigned long>;
  using VectorType = typename PointType::VectorType;

  /** Set the resolution level to be used for generating cells in the
   * Sphere. High values of this parameter will produce sphere with more
   * triangles. */
  /** Set/Get the name of the file to be read. */
  itkSetStringMacro(FileName);
  itkGetStringMacro(FileName);

protected:
  QuadEdgeMeshVTKPolyDataReader();
  ~QuadEdgeMeshVTKPolyDataReader() {}

  void
  PrintSelf(std::ostream & os, Indent indent) const override;

  /** Reads the file */
  void
  GenerateData() override;

  /** Filename to read */
  std::string m_FileName;
};
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#  include "itkQuadEdgeMeshVTKPolyDataReader.hxx"
#endif

#endif // _itkQuadEdgeMeshVTKPolyDataReader_h
