/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkPiecewiseRescaleQuadEdgeMeshFilter.h,v $
  Language:  C++
  Date:      $Date: 2011-05-03 13:09:05 $
  Version:   $Revision: 1.0 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkPiecewiseRescaleQuadEdgeMeshFilter_h
#define __itkPiecewiseRescaleQuadEdgeMeshFilter_h

#include "itkQuadEdgeMeshToQuadEdgeMeshFilter.h"

namespace itk
{
/**
 * \class PiecewiseRescaleQuadEdgeMeshFilter
 * \brief This filter rescales the scalar values of a mesh.
 *
 * The output mesh will have the same geometry as the input.
 *
 * It maps the input value that is in [min_in,cValue) to [min_out,cValue);
 *                                    (cValue,max_in] to (cValue,max_out].
 * and keeps cValue unchanged.
 *
 * \ingroup MeshFilters
 *
 */
template <class TInputMesh, class TOutputMesh>
class PiecewiseRescaleQuadEdgeMeshFilter :
  public QuadEdgeMeshToQuadEdgeMeshFilter<TInputMesh, TOutputMesh>
{
public:
  typedef PiecewiseRescaleQuadEdgeMeshFilter                        Self;
  typedef QuadEdgeMeshToQuadEdgeMeshFilter<TInputMesh, TOutputMesh> Superclass;
  typedef SmartPointer<Self>                                        Pointer;
  typedef SmartPointer<const Self>                                  ConstPointer;

  /** Run-time type information (and related methods).   */
  itkTypeMacro( PiecewiseRescaleQuadEdgeMeshFilter, QuadEdgeMeshToQuadEdgeMeshFilter );

  /** New macro for creation of through a Smart Pointer   */
  itkNewMacro( Self );

  typedef TInputMesh                                             InputMeshType;
  typedef typename InputMeshType::PixelType                      InputPixelType;
  typedef typename InputMeshType::PointDataContainer             InputPointDataContainer;
  typedef typename InputMeshType::PointDataContainerConstPointer InputPointDataContainerConstPointer;

  typedef TOutputMesh                                        OutputMeshType;
  typedef typename OutputMeshType::PixelType                 OutputPixelType;
  typedef typename OutputMeshType::PointDataContainer        OutputPointDataContainer;
  typedef typename OutputMeshType::PointDataContainerPointer OutputPointDataContainerPointer;

  /** Set/Get the mesh that will be deformed. */
  void SetInputMesh( const InputMeshType * mesh );

  const InputMeshType * GetInputMesh( void ) const;

  /** Set/Get min value of the output mesh scalars. */
  itkSetMacro( OutputMinimum, OutputPixelType );
  itkGetMacro( OutputMinimum, OutputPixelType );

  /** Set/Get max value of the output mesh scalars. */
  itkSetMacro( OutputMaximum, OutputPixelType );
  itkGetMacro( OutputMaximum, OutputPixelType );

  /** Get min value of the input mesh scalars. */
  itkGetMacro( InputMinimum, InputPixelType );

  /** Get max value of the input mesh scalars. */
  itkGetMacro( InputMaximum, InputPixelType );

  /** Set/Get value that we want to keep unchanged. */
  itkSetMacro( cValue, InputPixelType );
  itkGetMacro( cValue, InputPixelType );

  /** Get scale to transform mesh scalars in [min_in,cValue). */
  itkGetMacro( Scale_a, double );

  /** Get scale to transform mesh scalars in (cValue,max_in]. */
  itkGetMacro( Scale_b, double );
protected:
  PiecewiseRescaleQuadEdgeMeshFilter();
  ~PiecewiseRescaleQuadEdgeMeshFilter();

  void GenerateData();

private:

  PiecewiseRescaleQuadEdgeMeshFilter( const Self & ); // purposely not implemented
  void operator=( const Self & );                     // purposely not implemented

  OutputPixelType m_OutputMinimum;
  OutputPixelType m_OutputMaximum;

  InputPixelType m_InputMinimum;
  InputPixelType m_InputMaximum;

  InputPixelType m_cValue;

  double m_Scale_a;
  double m_Scale_b;
};
}

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkPiecewiseRescaleQuadEdgeMeshFilter.txx"
#endif

#endif
