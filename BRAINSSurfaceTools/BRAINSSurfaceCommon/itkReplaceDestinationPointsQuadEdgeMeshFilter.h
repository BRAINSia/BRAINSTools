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
#ifndef __itkReplaceDestinationPointsQuadEdgeMeshFilter_h
#define __itkReplaceDestinationPointsQuadEdgeMeshFilter_h

#include "itkQuadEdgeMeshToQuadEdgeMeshFilter.h"

namespace itk
{
/**
 * \class ReplaceDestinationPointsQuadEdgeMeshFilter
 * \brief This filter generates replace the coordinates of the input point mesh
 * with the coordinates of the input point set.
 *
 * This filter takes as input a PointSet, and a fixed Mesh, and assumes that
 * the points in the PointSet are one-to-one destination points for the points
 * in the fixed Mesh. Then, it replaces the coordinates of input mesh points
 * with the coordinates of input point set points.
 * \ingroup MeshFilters
 *
 */
template <class TInputMesh, class TInputPointSet>
class ReplaceDestinationPointsQuadEdgeMeshFilter :
  public QuadEdgeMeshToQuadEdgeMeshFilter<TInputMesh, TInputMesh>
{
public:
  typedef ReplaceDestinationPointsQuadEdgeMeshFilter Self;
  typedef QuadEdgeMeshToQuadEdgeMeshFilter<
      TInputMesh, TInputMesh>                            Superclass;
  typedef SmartPointer<Self>       Pointer;
  typedef SmartPointer<const Self> ConstPointer;

  /** Run-time type information (and related methods).   */
  itkTypeMacro( ReplaceDestinationPointsQuadEdgeMeshFilter, QuadEdgeMeshToQuadEdgeMeshFilter );

  /** New macro for creation of through a Smart Pointer   */
  itkNewMacro( Self );

  typedef TInputPointSet InputPointSetType;

  typedef TInputMesh InputMeshType;

  typedef typename Superclass::OutputMeshType  OutputMeshType;
  typedef typename Superclass::OutputPointType OutputPointType;

  /** Set Mesh whose grid defines the geometry and topology of the input PointSet.
   *  In a multi-resolution registration scenario, this will typically be the Input
   *  mesh at the current higher resolution level. */
  void SetInputMesh( const InputMeshType * inputMesh );

  const InputMeshType * GetInputMesh( void ) const;

  /** Set Mesh whose grid defines the geometry and topology of the input PointSet.
   *  In a multi-resolution registration scenario, this will typically be the Input
   *  mesh at the current higher resolution level. */
  void SetDestinationPoints( const InputPointSetType * destinationPointSet );

  const InputPointSetType * GetDestinationPoints( void ) const;

protected:
  ReplaceDestinationPointsQuadEdgeMeshFilter();
  ~ReplaceDestinationPointsQuadEdgeMeshFilter();

  void GenerateData();

private:

  ReplaceDestinationPointsQuadEdgeMeshFilter( const Self & ); // purposely not implemented
  void operator=( const Self & );                             // purposely not implemented
};
}

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkReplaceDestinationPointsQuadEdgeMeshFilter.hxx"
#endif

#endif
