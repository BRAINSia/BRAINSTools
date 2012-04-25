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
#ifndef __itkDeformationFieldFromTransformMeshFilter_h
#define __itkDeformationFieldFromTransformMeshFilter_h

#include "itkMeshToMeshFilter.h"
#include "itkTransform.h"

namespace itk
{
/**
 * \class DeformationFieldFromTransformMeshFilter
 * \brief Generate destination points from a Mesh and a Transform.
 *
 * This filter takes as input a Mesh and a Transform and produces as
 * output a point container with the list of Mesh points image resulting from
 * the Transform mapping.
 *
 * \ingroup MeshFilters
 *
 */
template <class TInputMesh, class TOutputMesh>
class DeformationFieldFromTransformMeshFilter :
  public MeshToMeshFilter<TInputMesh, TOutputMesh>
{
public:
  typedef DeformationFieldFromTransformMeshFilter Self;
  typedef MeshToMeshFilter<
      TInputMesh, TOutputMesh>                         Superclass;
  typedef SmartPointer<Self>       Pointer;
  typedef SmartPointer<const Self> ConstPointer;

  /** Run-time type information (and related methods).   */
  itkTypeMacro( DeformationFieldFromTransformMeshFilter, MeshToMeshFilter );

  /** New macro for creation of through a Smart Pointer   */
  itkNewMacro( Self );

  typedef TInputMesh                                           InputMeshType;
  typedef typename InputMeshType::ConstPointer                 InputMeshConstPointer;
  typedef typename InputMeshType::PointsContainerConstIterator InputPointsContainerConstIterator;

  typedef TOutputMesh                      OutputMeshType;
  typedef typename OutputMeshType::Pointer OutputMeshPointer;

  typedef typename OutputMeshType::PointsContainerIterator      OutputPointsContainerIterator;
  typedef typename OutputMeshType::PointsContainerConstIterator OutputPointsContainerConstIterator;

  itkStaticConstMacro( PointDimension, unsigned int, OutputMeshType::PointDimension );

  /** Transform typedef. */
  typedef Transform<double,
                    itkGetStaticConstMacro(PointDimension),
                    itkGetStaticConstMacro(PointDimension)>         TransformType;
  typedef typename TransformType::ConstPointer TransformPointerType;

  /** Set the coordinate transformation.  Set the coordinate transform that
   * will map the points of the input mesh to points of the output PointSet.
   * The points of the output PointSet are one-to-one the result of taking
   * points from the input Mesh and mapping them through the Transform.
   */
  itkSetConstObjectMacro( Transform, TransformType );

  /** Get a pointer to the coordinate transform. */
  itkGetConstObjectMacro( Transform, TransformType );
protected:
  DeformationFieldFromTransformMeshFilter();
  ~DeformationFieldFromTransformMeshFilter();

  void PrintSelf(std::ostream& os, Indent indent) const;

  void GenerateOutputInformation();

  void GenerateData();

private:

  DeformationFieldFromTransformMeshFilter( const Self & ); // purposely not implemented
  void operator=( const Self & );                          // purposely not implemented

  TransformPointerType m_Transform;             // Coordinate transform to use
};
}

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkDeformationFieldFromTransformMeshFilter.txx"
#endif

#endif
