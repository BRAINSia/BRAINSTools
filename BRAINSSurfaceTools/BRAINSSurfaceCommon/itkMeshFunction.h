/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkMeshFunction.h,v $
  Language:  C++
  Date:      $Date: 2008-10-17 13:35:26 $
  Version:   $Revision: 1.47 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkMeshFunction_h
#define __itkMeshFunction_h

#include "itkFunctionBase.h"
#include "itkMesh.h"
#include "itkPoint.h"

namespace itk
{
/** \class MeshFunction
 * \brief Evaluates a function of a Mesh at specified position.
 *
 * MeshFunction is a baseclass for all objects that evaluates
 * a function of a mesh at a given point.
 *
 * This class is templated over the input mesh type, the type
 * of the function output and the coordinate representation type
 * (e.g. float or double).
 *
 * The input mesh is set via method SetInputMesh().
 * Methods Evaluate, EvaluateAtIndex and EvaluateAtContinuousIndex
 * respectively evaluates the function at an geometric point,
 * mesh index and continuous mesh index.
 *
 * \warning Mesh information is cached during SetInputMesh( mesh ). If the mesh
 * has changed one must call SetInputMesh( mesh ) again to update the cache to
 * the current values. The cached values are the ones of the KdTree point locator.
 *
 * \sa Point
 *
 * \ingroup MeshFunctions
 */
template <class TInputMesh, class TOutput>
class ITK_EXPORT MeshFunction :
  public         FunctionBase<typename TInputMesh::PointType, TOutput>
{
public:
  /** Dimension underlying input mesh. */
  itkStaticConstMacro(MeshDimension, unsigned int, TInputMesh::PointDimension);

  /** Standard class typedefs. */
  typedef MeshFunction Self;
  typedef FunctionBase<
      typename TInputMesh::PointType, TOutput>           Superclass;
  typedef SmartPointer<Self>       Pointer;
  typedef SmartPointer<const Self> ConstPointer;

  /** Run-time type information (and related methods). */
  itkTypeMacro(MeshFunction, FunctionBase);

  /** InputMeshType typedef support. */
  typedef TInputMesh InputMeshType;

  /** InputPixel typedef support */
  typedef typename InputMeshType::PixelType InputPixelType;

  /** InputMeshPointer typedef support */
  typedef typename InputMeshType::ConstPointer InputMeshConstPointer;

  /** OutputType typedef support. */
  typedef TOutput OutputType;

  /** CoordRepType typedef support. */
  typedef typename InputMeshType::CoordRepType CoordRepType;

  /** Point Type. */
  typedef typename InputMeshType::PointType PointType;

  /** Set the input mesh.
   * \warning this method caches information.  If the Mesh has changed, user
   * must call SetInputMesh again to update cached values. */
  virtual void SetInputMesh( const InputMeshType * ptr );

  /** Get the input mesh. */
  const InputMeshType * GetInputMesh() const;

  /** Evaluate the function at specified Point position.
   * Subclasses must provide this method. */
  virtual TOutput Evaluate( const PointType& point ) const = 0;

protected:
  MeshFunction();
  ~MeshFunction()
  {
  }

  void PrintSelf(std::ostream& os, Indent indent) const;

  /** Const pointer to the input image. */
  InputMeshConstPointer m_Mesh;
private:
  MeshFunction(const Self &);   // purposely not implemented
  void operator=(const Self &); // purposely not implemented
};
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkMeshFunction.hxx"
#endif

#endif
