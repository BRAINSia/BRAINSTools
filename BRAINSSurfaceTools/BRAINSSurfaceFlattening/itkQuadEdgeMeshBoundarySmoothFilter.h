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
#ifndef __itkQuadEdgeMeshBoundarySmoothFilter_h
#define __itkQuadEdgeMeshBoundarySmoothFilter_h

#include <itkQuadEdgeMeshToQuadEdgeMeshFilter.h>

namespace itk
{
/**
 * \class QuadEdgeMeshBoundarySmoothFilter
 * \brief This filter takes as inout two surfaces that represent
 * a single genus 0 surface. The smoothing will generate a smooth
 * boundary between the two surfaces. The filter produces two outputs
 * that are the resulting smoothed surfaces.
 *
 * The filter is supposed to be used after itkQuadEdgeMeshSplitFilter
 * and before mapping the surface to sphere.

 * itkQuadEdgeMeshSplitFilter takes one surface as input and
 * split it into two surfaces. There are "teeth" left on boundaries
 * of splitted surfaces because of the triangulars on the surface.
 *
 * The smoothing method is repeated for the user specified number
 * of iterations or until there is no triangles to move between
 * surfaces.
 *
 * Note: the total number of triangles on both inputs have to be
 * the same after each iteration. Missing cell or extra cell to
 * the whole surface is not allowed.
 *
 * The input mesh must share the same mesh type, but the output
 * may differ.
 * \ingroup MeshFilters
 *
 * \sa itkQuadEdgeMeshSplitFilter
 *
 */
template <class TInputMesh, class TOutputMesh>
class QuadEdgeMeshBoundarySmoothFilter :
  public QuadEdgeMeshToQuadEdgeMeshFilter<TInputMesh, TOutputMesh>
{
public:
  typedef QuadEdgeMeshBoundarySmoothFilter                          Self;
  typedef SmartPointer<Self>                                        Pointer;
  typedef SmartPointer<const Self>                                  ConstPointer;
  typedef QuadEdgeMeshToQuadEdgeMeshFilter<TInputMesh, TOutputMesh> Superclass;

  /** Input types. */
  typedef TInputMesh                              InputMeshType;
  typedef typename InputMeshType::Pointer         InputMeshPointer;
  typedef typename InputMeshType::ConstPointer    InputMeshConstPointer;
  typedef typename InputMeshType::CoordRepType    InputCoordRepType;
  typedef typename InputMeshType::PointType       InputPointType;
  typedef typename InputMeshType::PointIdentifier InputPointIdentifier;
  typedef typename InputMeshType::QEType          InputQEType;
  typedef typename InputMeshType::QEPrimal        InputQEPrimal;
  typedef typename InputMeshType::VectorType      InputVectorType;
  typedef typename InputMeshType::CellType        InputCellType;
  typedef typename InputMeshType::CellIdentifier  InputCellIdentifier;

  typedef typename InputMeshType::PointsContainerConstPointer
    InputPointsContainerConstPointer;
  typedef typename InputMeshType::PointsContainerConstIterator
    InputPointsContainerConstIterator;
  typedef typename InputMeshType::CellsContainerConstPointer
    InputCellsContainerConstPointer;
  typedef typename InputMeshType::CellsContainerConstIterator
    InputCellsContainerConstIterator;

  typedef typename InputMeshType::EdgeCellType    InputEdgeCellType;
  typedef typename InputMeshType::PolygonCellType InputPolygonCellType;
  typedef typename InputMeshType::PointIdList     InputPointIdList;
  typedef typename InputMeshType::CellTraits      InputCellTraits;
  typedef typename InputCellTraits::PointIdInternalIterator
    InputPointsIdInternalIterator;
//    typedef typename InputQEPrimal::IteratorGeom    InputQEIterator;

  /** Output types. */
  typedef TOutputMesh                              OutputMeshType;
  typedef typename OutputMeshType::Pointer         OutputMeshPointer;
  typedef typename OutputMeshType::ConstPointer    OutputMeshConstPointer;
  typedef typename OutputMeshType::CoordRepType    OutputCoordRepType;
  typedef typename OutputMeshType::PointType       OutputPointType;
  typedef typename OutputMeshType::PointIdentifier OutputPointIdentifier;
  typedef typename OutputMeshType::PointIdList     OutputPointIdListType;

  typedef typename OutputMeshType::EdgeListPointerType EdgeListPointerType;
  typedef typename OutputMeshType::QEPrimal            OutputQEPrimal;
  typedef typename OutputQEPrimal::IteratorGeom        OutputQEIterator;
  typedef typename OutputMeshType::PolygonCellType     OutputPolygonCellType;
  typedef typename OutputMeshType::PointsContainerPointer
    OutputPointsContainerPointer;
  typedef typename OutputMeshType::PointsContainerIterator
    OutputPointsContainerIterator;
  typedef typename OutputMeshType::CellsContainerConstPointer
    OutputCellsContainerConstPointer;
  typedef typename OutputMeshType::CellsContainerConstIterator
    OutputCellsContainerConstIterator;
public:
  itkNewMacro( Self );
  itkTypeMacro( QuadEdgeMeshBoundarySmoothFilter, QuadEdgeMeshToQuadEdgeMeshFilter );

  /** Set/Get the first input mesh */
  void SetInputMesh1( const InputMeshType * mesh1 );

  const InputMeshType * GetInputMesh1( void ) const;

  /** Set/Get the second input mesh */
  void SetInputMesh2( const InputMeshType * mesh2 );

  const InputMeshType * GetInputMesh2( void ) const;

  /** Get the first smoothed hemisphere */
  OutputMeshType * GetOutputMesh1( void );

  /** Get the second smoothed hemisphere */
  OutputMeshType * GetOutputMesh2( void );

  /** Set/Get the number of iterations. */
  itkSetMacro( Iterations, int );
  itkGetMacro( Iterations, int );
protected:
  QuadEdgeMeshBoundarySmoothFilter();
  ~QuadEdgeMeshBoundarySmoothFilter();

  void CopyInputMeshesToOutputMeshes();

  int AdjustBoundary( OutputMeshType * deleteMesh, OutputMeshType * addMesh);

  virtual void GenerateData() override;

private:
  QuadEdgeMeshBoundarySmoothFilter( const Self & );
  void operator =( const Self & );

  int m_Iterations;
};
}

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkQuadEdgeMeshBoundarySmoothFilter.hxx"
#endif

#endif
