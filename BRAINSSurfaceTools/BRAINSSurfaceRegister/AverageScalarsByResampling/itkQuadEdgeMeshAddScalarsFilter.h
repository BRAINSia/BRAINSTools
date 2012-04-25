/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkQuadEdgeMeshAddScalarsFilter.h,v $
  Language:  C++
  Date:      $Date: 2008-05-14 09:26:05 $
  Version:   $Revision: 1.21 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkQuadEdgeMeshAddScalarsFilter_h
#define __itkQuadEdgeMeshAddScalarsFilter_h

#include "itkQuadEdgeMeshToQuadEdgeMeshFilter.h"

namespace itk
{
/**
 * \class QuadEdgeMeshAddScalarsFilter
 * \brief add scalars from two input meshes together
 *
 * This filter takes inputMesh1 and inputMesh2
 * as inputs and generates the outputMesh that has
 * scalars as the sum of that of inputMesh1 and inputMesh2.
 *
 * It requires that input meshes have the point
 * -point correspondence.
 *
 * \ingroup MeshFilters
 *
 */
template <class TInputMesh1, class TInputMesh2, class TOutputMesh>
class QuadEdgeMeshAddScalarsFilter :
  public QuadEdgeMeshToQuadEdgeMeshFilter<TInputMesh1, TInputMesh1>
{
public:
  typedef QuadEdgeMeshAddScalarsFilter Self;
  typedef QuadEdgeMeshToQuadEdgeMeshFilter<
      TInputMesh1, TInputMesh1>                         Superclass;
  typedef SmartPointer<Self>       Pointer;
  typedef SmartPointer<const Self> ConstPointer;

  /** Run-time type information (and related methods).   */
  itkTypeMacro( QuadEdgeMeshAddScalarsFilter, QuadEdgeMeshToQuadEdgeMeshFilter );

  /** New macro for creation of through a Smart Pointer   */
  itkNewMacro( Self );

  typedef TInputMesh1 InputMeshType1;

  typedef TInputMesh2 InputMeshType2;

  typedef TOutputMesh OutputMeshType;

  /** Set/Get the input mesh 1. */
  void SetInput1( const InputMeshType1 * mesh );

  const InputMeshType1 * GetInput1( void ) const;

  /** Set/Get the input mesh2. */
  void SetInput2( const InputMeshType2 * mesh );

  const InputMeshType2 * GetInput2( void ) const;

protected:
  QuadEdgeMeshAddScalarsFilter();
  ~QuadEdgeMeshAddScalarsFilter();

  void GenerateData();

private:

  QuadEdgeMeshAddScalarsFilter( const Self & ); // purposely not implemented
  void operator=( const Self & );               // purposely not implemented
};
}

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkQuadEdgeMeshAddScalarsFilter.txx"
#endif

#endif
