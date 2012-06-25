/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkConvertVTKToQuadEdgeMeshFilter.h,v $
  Language:  C++
  Date:      $Date: 2008-01-15 19:10:40 $
  Version:   $Revision: 1.8 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
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
  void GenerateData();

private:
  ConvertVTKToQuadEdgeMeshFilter(const Self &); // purposely not implemented
  void operator=(const Self &);                 // purposely not implemented

  vtkPolyData * m_inputPolyData;
};
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkConvertVTKToQuadEdgeMeshFilter.hxx"
#endif

#endif // _itkConvertVTKToQuadEdgeMeshFilter_h
