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
#ifndef __itkConvertVTKToQuadEdgeMeshFilter_h
#define __itkConvertVTKToQuadEdgeMeshFilter_h

#include "itkMesh.h"
#include "itkMeshSource.h"
#include "itkTriangleCell.h"
#include "itkMapContainer.h"

#include "vtkPolyData.h"
#include "vtkPointData.h"
#include "vtkIdList.h"

namespace itk
{
/** \class ConvertVTKToQuadEdgeMeshFilter
 * \brief
 * take a vtkPolyData as input and create an itkMesh.
 *
 * Caveat: itkConvertVTKToQuadEdgeMeshFilter can only convert triangle meshes.
 *         Use vtkTriangleFilter to convert your mesh to a triangle mesh.
 */
template <class TOutputMesh>
class ConvertVTKToQuadEdgeMeshFilter : public MeshSource<TOutputMesh>
{
public:
  /** Standard "Self" typedef. */
  typedef ConvertVTKToQuadEdgeMeshFilter Self;
  typedef MeshSource<TOutputMesh>        Superclass;
  typedef SmartPointer<Self>             Pointer;
  typedef SmartPointer<const Self>       ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(ConvertVTKToQuadEdgeMeshFilter, MeshSource);

  /** Hold on to the type information specified by the template parameters. */
  typedef TOutputMesh                        OutputMeshType;
  typedef typename OutputMeshType::PointType PointType;
  typedef typename OutputMeshType::PixelType PixelType;

  /** Some convenient typedefs. */
  typedef typename OutputMeshType::Pointer         OutputMeshPointer;
  typedef typename OutputMeshType::CellTraits      CellTraits;
  typedef typename OutputMeshType::CellIdentifier  CellIdentifier;
  typedef typename OutputMeshType::CellType        CellType;
  typedef typename OutputMeshType::CellAutoPointer CellAutoPointer;
  typedef typename OutputMeshType::PointIdentifier PointIdentifier;
  typedef typename OutputMeshType::PointIdList     PointIdListType;
  typedef typename CellTraits::PointIdIterator     PointIdIterator;

  /** Define the triangular cell types which form the surface  */
  typedef TriangleCell<CellType> TriangleCellType;

  typedef typename TriangleCellType::SelfAutoPointer
    TriangleCellAutoPointer;

  typedef std::pair<unsigned long, unsigned long>    IndexPairType;
  typedef MapContainer<IndexPairType, unsigned long> PointMapType;
  typedef typename PointType::VectorType             VectorType;

  /* set/get the input polydata */
  void SetPolyData( vtkPolyData * );

  vtkPolyData * GetPolyData();

protected:
  ConvertVTKToQuadEdgeMeshFilter();
  ~ConvertVTKToQuadEdgeMeshFilter()
  {
  };

  /** convert the polydata */
  void GenerateData() ITK_OVERRIDE;

private:
  ITK_DISALLOW_COPY_AND_ASSIGN(ConvertVTKToQuadEdgeMeshFilter);

  vtkPolyData * m_inputPolyData;
};
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkConvertVTKToQuadEdgeMeshFilter.hxx"
#endif

#endif // _itkConvertVTKToQuadEdgeMeshFilter_h
