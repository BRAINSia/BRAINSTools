/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkNodeScalarGradientCalculator.h,v $
  Language:  C++
  Date:      $Date: 2008-05-14 09:26:05 $
  Version:   $Revision: 1.21 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkAssignScalarValuesQuadEdgeMeshFilter_h
#define __itkAssignScalarValuesQuadEdgeMeshFilter_h

#include "itkQuadEdgeMeshToQuadEdgeMeshFilter.h"

namespace itk
{
/**
 * \class AssignScalarValuesQuadEdgeMeshFilter
 * \brief This filter was initially proposed to assign resampled moving values back onto FixedMesh.
 *
 *
 * This filter takes two inputs as source and input meshes and assumes that
 * the points in the input are one-to-one for the points
 * in the source Mesh. Then, it assigns the scalar values of the source mesh points
 * to the points of the input mesh.
 * \ingroup MeshFilters
 *
 */
template <class TInputMesh, class TSourceMesh, class TOutputMesh>
class AssignScalarValuesQuadEdgeMeshFilter :
  public QuadEdgeMeshToQuadEdgeMeshFilter<TInputMesh, TOutputMesh>
{
public:
  typedef AssignScalarValuesQuadEdgeMeshFilter Self;
  typedef QuadEdgeMeshToQuadEdgeMeshFilter<
      TInputMesh, TOutputMesh>                     Superclass;
  typedef SmartPointer<Self>       Pointer;
  typedef SmartPointer<const Self> ConstPointer;

  /** Run-time type information (and related methods).   */
  itkTypeMacro( AssignScalarValuesQuadEdgeMeshFilter, QuadEdgeMeshToQuadEdgeMeshFilter );

  /** New macro for creation of through a Smart Pointer   */
  itkNewMacro( Self );

  typedef TSourceMesh                                  SourceMeshType;
  typedef typename  SourceMeshType::Pointer            SourceMeshPointer;
  typedef typename  SourceMeshType::PointDataContainer SourcePointDataContainer;

  typedef typename  Superclass::InputMeshType InputMeshType;

  typedef typename  Superclass::OutputMeshType           OutputMeshType;
  typedef typename  OutputMeshType::Pointer              OutputMeshPointer;
  typedef typename  Superclass::OutputPointDataContainer OutputPointDataContainer;
  typedef typename  OutputPointDataContainer::Pointer    OutputPointDataContainerPointer;
  typedef typename  OutputPointDataContainer::Iterator   OutputPointDataContainerIterator;

  // Set source mesh who has the scalar values we want.
  void SetSourceMesh( const SourceMeshType * sourceMesh );

  const SourceMeshType * GetSourceMesh( void ) const;

  // Set input mesh who wants to have scalar values from sourceMesh.
  void SetInputMesh( const InputMeshType * inputMesh );

  const InputMeshType * GetInputMesh( void ) const;

protected:
  AssignScalarValuesQuadEdgeMeshFilter();
  ~AssignScalarValuesQuadEdgeMeshFilter();

  void GenerateData();

private:

  AssignScalarValuesQuadEdgeMeshFilter( const Self & ); // purposely not implemented
  void operator=( const Self & );                       // purposely not implemented
};
}

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkAssignScalarValuesQuadEdgeMeshFilter.txx"
#endif

#endif
