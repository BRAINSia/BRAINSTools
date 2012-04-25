/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkScalarQuadEdgeMeshToListAdaptor.h,v $
  Language:  C++
  Date:      $Date: 2008-05-14 09:26:05 $
  Version:   $Revision: 1.21 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkScalarQuadEdgeMeshToListAdaptor_h
#define __itkScalarQuadEdgeMeshToListAdaptor_h

#include "itkQuadEdgeMeshToQuadEdgeMeshFilter.h"

#include "itkFixedArray.h"
#include "itkListSample.h"

namespace itk
{
/**
 * \class ScalarQuadEdgeMeshToListAdaptor
 * \brief The adaptor takes one input mesh and generate
 *  a ListSample of scalar values
 *
 *
 *
 * \ingroup MeshFilters
 *
 */
template <class TInputMesh>
class ScalarQuadEdgeMeshToListAdaptor :
  public QuadEdgeMeshToQuadEdgeMeshFilter<TInputMesh, TInputMesh>
{
public:
  typedef ScalarQuadEdgeMeshToListAdaptor                          Self;
  typedef QuadEdgeMeshToQuadEdgeMeshFilter<TInputMesh, TInputMesh> Superclass;
  typedef SmartPointer<Self>                                       Pointer;
  typedef SmartPointer<const Self>                                 ConstPointer;

  /** Run-time type information (and related methods).   */
  itkTypeMacro( ScalarQuadEdgeMeshToListAdaptor, QuadEdgeMeshToQuadEdgeMeshFilter );

  /** New macro for creation of through a Smart Pointer   */
  itkNewMacro( Self );

  typedef TInputMesh                                        InputMeshType;
  typedef typename InputMeshType::PointType                 InputPointType;
  typedef typename InputMeshType::PixelType                 InputPixelType;
  typedef typename InputMeshType::PointIdentifier           PointIdentifier;
  typedef typename InputMeshType::PointDataContainer        InputPointDataContainer;
  typedef typename InputMeshType::PointDataContainerPointer InputPointDataContainerPointer;

  typedef typename InputMeshType::PixelType                           MeasurementType;
  typedef typename itk::FixedArray<MeasurementType, 1>                MeasurementVectorType;
  typedef typename itk::Statistics::ListSample<MeasurementVectorType> ListSampleType;

  typedef typename ListSampleType::Pointer ListSamplePointerType;

  /** Set/Get the mesh. */
  void SetMesh( const InputMeshType * mesh );

  const InputMeshType * GetMesh( void ) const;

  /** Get the calculated ListSample. */
  itkGetObjectMacro( Sample, ListSampleType );

  void Compute();

protected:
  ScalarQuadEdgeMeshToListAdaptor();
  ~ScalarQuadEdgeMeshToListAdaptor();
private:

  ScalarQuadEdgeMeshToListAdaptor( const Self & ); // purposely not implemented
  void operator=( const Self & );                  // purposely not implemented

  ListSamplePointerType m_Sample;
};
}

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkScalarQuadEdgeMeshToListAdaptor.txx"
#endif

#endif
