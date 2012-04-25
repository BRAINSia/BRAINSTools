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
#ifndef __itkWarpQuadEdgeMeshFilter_h
#define __itkWarpQuadEdgeMeshFilter_h

#include "itkQuadEdgeMeshToQuadEdgeMeshFilter.h"

#include "itkInterpolateMeshFunction.h"
#include "itkLinearInterpolateMeshFunction.h"

namespace itk
{
/**
 * \class WarpQuadEdgeMeshFilter
 * \brief This filter warps the mesh of its first input by applying the
 * deformation field and scalar values from the reference mesh.
 *
 * This filter will take three inputs.
 * Input 0: input fixed mesh.
 * Input 1: reference moving mesh.
 * Input 2: deformation field.
 * The deformation field and the input fixed mesh
 * need to have one to one correspondence.
 *
 * The output mesh will be the copy of the input fixed mesh.
 * Scalar values of the Points of the output mesh will be got
 * by interpolating the deformed points on reference mesh after
 * applying the displacement vector from deformation field.
 *
 * \ingroup MeshFilters
 *
 */
template <class TInputMesh, class TReferenceMesh, class TDeformationField>
class WarpQuadEdgeMeshFilter :
  public QuadEdgeMeshToQuadEdgeMeshFilter<TInputMesh, TInputMesh>
{
public:
  typedef WarpQuadEdgeMeshFilter                                   Self;
  typedef QuadEdgeMeshToQuadEdgeMeshFilter<TInputMesh, TInputMesh> Superclass;
  typedef SmartPointer<Self>                                       Pointer;
  typedef SmartPointer<const Self>                                 ConstPointer;

  /** Run-time type information (and related methods).   */
  itkTypeMacro( WarpQuadEdgeMeshFilter, QuadEdgeMeshToQuadEdgeMeshFilter );

  /** New macro for creation of through a Smart Pointer   */
  itkNewMacro( Self );

  typedef TInputMesh InputMeshType;

  typedef TReferenceMesh ReferenceMeshType;

  typedef TInputMesh                                         OutputMeshType;
  typedef typename OutputMeshType::PointType                 OutputPointType;
  typedef typename OutputMeshType::PointsContainer           OutputPointsContainer;
  typedef typename OutputMeshType::PointsContainerPointer    OutputPointsContainerPointer;
  typedef typename OutputMeshType::PointDataContainer        OutputPointDataContainer;
  typedef typename OutputMeshType::PointDataContainerPointer OutputPointDataContainerPointer;

  typedef TDeformationField                           DeformationFieldType;
  typedef typename DeformationFieldType::ConstPointer DeformationFieldPointer;
  typedef typename DeformationFieldType::PixelType    DisplacementType;

  typedef typename DeformationFieldType::PointDataContainer        DisplacementVectorContainer;
  typedef typename DeformationFieldType::PointDataContainerPointer DisplacementVectorContainerPointer;

  /** Interpolator typedef. */
  typedef InterpolateMeshFunction<InputMeshType>       InterpolatorType;
  typedef typename InterpolatorType::Pointer           InterpolatorPointerType;
  typedef LinearInterpolateMeshFunction<InputMeshType> DefaultInterpolatorType;

  /** Set/Get the mesh that will be deformed. */
  void SetInputMesh( const InputMeshType * mesh );

  const InputMeshType * GetInputMesh( void ) const;

  /** Set/Get the mesh that carried the deformation field as pixel data. */
  void SetReferenceMesh( const ReferenceMeshType * mesh );

  const ReferenceMeshType * GetReferenceMesh( void ) const;

  /** Set/Get the mesh that carried the deformation field as pixel data. */
  void SetDeformationField( const DeformationFieldType * field );

  const DeformationFieldType * GetDeformationField( void ) const;

  /** Set the interpolator function.  The default is a linear interpolator. */
  itkSetObjectMacro( Interpolator, InterpolatorType );

  /** Get a pointer to the interpolator function. */
  itkGetConstObjectMacro( Interpolator, InterpolatorType );
protected:
  WarpQuadEdgeMeshFilter();
  ~WarpQuadEdgeMeshFilter();

  void GenerateData();

private:

  WarpQuadEdgeMeshFilter( const Self & ); // purposely not implemented
  void operator=( const Self & );         // purposely not implemented

  InterpolatorPointerType m_Interpolator;
};
}

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkWarpQuadEdgeMeshFilter.txx"
#endif

#endif
