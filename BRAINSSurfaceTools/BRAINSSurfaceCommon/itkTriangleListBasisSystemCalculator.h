/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkTriangleListBasisSystemCalculator.h,v $
  Language:  C++
  Date:      $Date: 2008-10-17 13:35:26 $
  Version:   $Revision: 1.47 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkTriangleListBasisSystemCalculator_h
#define __itkTriangleListBasisSystemCalculator_h

#include "itkMesh.h"
#include "itkPoint.h"
#include "itkTriangleBasisSystemCalculator.h"

namespace itk
{
/** \class TriangleListBasisSystemCalculator
 * \brief  Computes basis coefficients for a list of triangular cells.
 *
 * TriangleListBasisSystemCalculator computes basis coefficients
 * within a list of triangles contained in a mesh. Basis coefficients
 * can be used thereafter for interpolation and gradient computation
 * within those triangles.
 *
 * This class is templated over the input vector type and dimension of basis.
 *
 * \sa TriangleBasisSystem
 *
 */
template <class TMesh, class TBasisSystem>
class TriangleListBasisSystemCalculator : public Object
{
public:
  /** Standard class typedefs. */
  typedef TriangleListBasisSystemCalculator Self;
  typedef Object                            Superclass;
  typedef SmartPointer<Self>                Pointer;
  typedef SmartPointer<const Self>          ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Standard part of every itk Object. */
  itkTypeMacro(TriangleListBasisSystemCalculator, Object);

  typedef TBasisSystem                                  BasisSystemType;
  typedef typename TMesh::CellIdentifier                CellIdentifier;
  typedef VectorContainer<CellIdentifier, TBasisSystem> BasisSystemListType;
  typedef typename BasisSystemListType::ConstIterator   BasisSystemListIterator;
  typedef typename BasisSystemListType::Pointer         BasisSystemListPointer;

  typedef TMesh                                  MeshType;
  typedef typename MeshType::PointType           PointType;
  typedef typename MeshType::CellType            CellType;
  typedef typename PointType::VectorType         VectorType;
  typedef typename MeshType::ConstPointer        MeshConstPointer;
  typedef typename MeshType::CellsContainer      CellsContainer;
  typedef typename CellsContainer::Iterator      CellsContainerIterator;
  typedef typename CellsContainer::ConstIterator CellsContainerConstIterator;

  /** Set/Get the input mesh. */
  itkSetConstObjectMacro( InputMesh, MeshType );
  itkGetConstObjectMacro( InputMesh, MeshType );

  /** Compute the basis system at every triangle. */
  void Calculate();

  /** Get the list of basis systems. */
  itkGetConstObjectMacro( BasisSystemList, BasisSystemListType );
protected:
  TriangleListBasisSystemCalculator();
  virtual ~TriangleListBasisSystemCalculator();
private:
  MeshConstPointer       m_InputMesh;
  BasisSystemListPointer m_BasisSystemList;
};
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkTriangleListBasisSystemCalculator.hxx"
#endif

#endif
