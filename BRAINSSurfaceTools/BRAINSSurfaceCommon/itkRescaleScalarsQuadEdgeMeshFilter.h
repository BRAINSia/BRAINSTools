/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkRescaleScalarsQuadEdgeMeshFilter.h,v $
  Language:  C++
  Date:      $Date: 2010-05-14 09:26:05 $
  Version:   $Revision: 1.21 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkRescaleScalarsQuadEdgeMeshFilter_h
#define __itkRescaleScalarsQuadEdgeMeshFilter_h

#include "itkQuadEdgeMeshToQuadEdgeMeshFilter.h"

namespace itk
{
/**
 * \class RescaleScalarsQuadEdgeMeshFilter
 * \brief This filter rescales the scalar values of a mesh.
 *
 * The output mesh will have the same geometry as the input.
 *
 *
 * It searches min and max values of the input mesh scalars automatically.
 * Applies a linear transformation to the scalar levels of the input Image.
 * by min and max values given by users.
 *
 * \ingroup MeshFilters
 *
 */
template <class TMesh>
class RescaleScalarsQuadEdgeMeshFilter :
  public QuadEdgeMeshToQuadEdgeMeshFilter<TMesh, TMesh>
{
public:
  typedef RescaleScalarsQuadEdgeMeshFilter               Self;
  typedef QuadEdgeMeshToQuadEdgeMeshFilter<TMesh, TMesh> Superclass;
  typedef SmartPointer<Self>                             Pointer;
  typedef SmartPointer<const Self>                       ConstPointer;

  /** Run-time type information (and related methods).   */
  itkTypeMacro( RescaleScalarsQuadEdgeMeshFilter, QuadEdgeMeshToQuadEdgeMeshFilter );

  /** New macro for creation of through a Smart Pointer   */
  itkNewMacro( Self );

  typedef TMesh                                                  InputMeshType;
  typedef typename InputMeshType::PixelType                      InputPixelType;
  typedef typename InputMeshType::PointDataContainer             InputPointDataContainer;
  typedef typename InputMeshType::PointDataContainerConstPointer InputPointDataContainerConstPointer;

  typedef TMesh                                              OutputMeshType;
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

  /** Get scale to transform mesh scalars. */
  itkGetMacro( Scale, double );
protected:
  RescaleScalarsQuadEdgeMeshFilter();
  ~RescaleScalarsQuadEdgeMeshFilter();

  void GenerateData();

private:

  RescaleScalarsQuadEdgeMeshFilter( const Self & ); // purposely not implemented
  void operator=( const Self & );                   // purposely not implemented

  OutputPixelType m_OutputMinimum;
  OutputPixelType m_OutputMaximum;

  InputPixelType m_InputMinimum;
  InputPixelType m_InputMaximum;

  double m_Scale;
};
}

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkRescaleScalarsQuadEdgeMeshFilter.hxx"
#endif

#endif
