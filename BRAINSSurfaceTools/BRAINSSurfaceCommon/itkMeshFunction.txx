/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkMeshFunction.txx,v $
  Language:  C++
  Date:      $Date: 2008-10-17 13:35:26 $
  Version:   $Revision: 1.17 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkMeshFunction_txx
#define __itkMeshFunction_txx

#include "itkMeshFunction.h"

namespace itk
{
/**
 * Constructor
 */
template <class TInputMesh, class TOutput>
MeshFunction<TInputMesh, TOutput>
::MeshFunction()
{
  m_Mesh = NULL;
}

/**
 * Standard "PrintSelf" method
 */
template <class TInputMesh, class TOutput>
void
MeshFunction<TInputMesh, TOutput>
::PrintSelf( std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf( os, indent );

  os << indent << "InputMesh: " << m_Mesh.GetPointer() << std::endl;
}

/**
 * Initialize by setting the input mesh
 */
template <class TInputMesh, class TOutput>
void
MeshFunction<TInputMesh, TOutput>
::SetInputMesh( const InputMeshType * ptr )
{
  this->m_Mesh = ptr;

  if( ptr )
    {
    // FIXME Add here the point locator...
    }
}

/**
 * Return the input mesh
 */
template <class TInputMesh, class TOutput>
const typename MeshFunction<TInputMesh, TOutput>::InputMeshType
* MeshFunction<TInputMesh, TOutput>
::GetInputMesh() const
  {
  return m_Mesh;
  }
} // end namespace itk

#endif
