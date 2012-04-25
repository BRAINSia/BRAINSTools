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
#ifndef __itkNormalizeScalarsQuadEdgeMeshFilter_h
#define __itkNormalizeScalarsQuadEdgeMeshFilter_h

#include "itkQuadEdgeMeshToQuadEdgeMeshFilter.h"

namespace itk
{
/**
 * \class NormalizeScalarsQuadEdgeMeshFilter
 * \brief This filter normalizes the scalar values of a mesh.
 *
 * The output mesh will have the same geometry.
 *
 * The normalization is done by first subtracting the median value and the
 * dividing by the standard deviation of the resulting values. This process is
 * applied multiple times under the control of a user-provided number of
 * iteraions.
 *
 * \ingroup MeshFilters
 *
 */
template <class TMesh>
class NormalizeScalarsQuadEdgeMeshFilter :
  public QuadEdgeMeshToQuadEdgeMeshFilter<TMesh, TMesh>
{
public:
  typedef NormalizeScalarsQuadEdgeMeshFilter             Self;
  typedef QuadEdgeMeshToQuadEdgeMeshFilter<TMesh, TMesh> Superclass;
  typedef SmartPointer<Self>                             Pointer;
  typedef SmartPointer<const Self>                       ConstPointer;

  /** Run-time type information (and related methods).   */
  itkTypeMacro( NormalizeScalarsQuadEdgeMeshFilter, QuadEdgeMeshToQuadEdgeMeshFilter );

  /** New macro for creation of through a Smart Pointer   */
  itkNewMacro( Self );

  typedef TMesh                                              InputMeshType;
  typedef TMesh                                              OutputMeshType;
  typedef typename OutputMeshType::PointDataContainer        OutputPointDataContainer;
  typedef typename OutputMeshType::PointDataContainerPointer OutputPointDataContainerPointer;

  /** Set/Get the mesh that will be deformed. */
  void SetInputMesh( const InputMeshType * mesh );

  const InputMeshType * GetInputMesh( void ) const;

  /** Set/Get number of iterations of normalization. */
  itkSetMacro( NumberOfIterations, unsigned int );
  itkGetMacro( NumberOfIterations, unsigned int );
protected:
  NormalizeScalarsQuadEdgeMeshFilter();
  ~NormalizeScalarsQuadEdgeMeshFilter();

  void GenerateData();

private:

  NormalizeScalarsQuadEdgeMeshFilter( const Self & ); // purposely not implemented
  void operator=( const Self & );                     // purposely not implemented

  unsigned int m_NumberOfIterations;
};
}

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkNormalizeScalarsQuadEdgeMeshFilter.txx"
#endif

#endif
