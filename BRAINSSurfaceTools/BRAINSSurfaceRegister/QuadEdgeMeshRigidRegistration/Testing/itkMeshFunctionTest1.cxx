/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkFreeSurferBinarySurfaceReaderQuadEdgeMeshTest.cxx,v $
  Language:  C++
  Date:      $Date: 2008-06-16 02:07:17 $
  Version:   $Revision: 1.2 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/

#ifdef _MSC_VER
#pragma warning ( disable : 4786 )
#endif

#include "itkQuadEdgeMesh.h"
#include "itkMeshFunction.h"

#include <iostream>

namespace itk
{
template <class TInputMesh, class TOutput>
class MeshFunctionTestHelper : public MeshFunction<TInputMesh, TOutput>
{
public:
  typedef MeshFunctionTestHelper            Self;
  typedef MeshFunction<TInputMesh, TOutput> Superclass;
  typedef SmartPointer<Self>                Pointer;
  typedef SmartPointer<const Self>          ConstPointer;

  itkTypeMacro(MeshFunctionTestHelper, MeshFunction);

  itkNewMacro( Self );

  typedef typename Superclass::OutputType OutputType;
  typedef typename Superclass::PointType  PointType;

  virtual OutputType Evaluate( const PointType & itkNotUsed(point) ) const
  {
    return OutputType();
  }

  void RunTest()
  {
    std::cout << "Superclass name = " << this->Superclass::GetNameOfClass() << std::endl;
    std::cout << "Class name      = " << this->GetNameOfClass() << std::endl;
  }
};
}

int main(int, char * [] )
{
  typedef itk::QuadEdgeMesh<float, 3> MeshType;

  MeshType::Pointer mesh = MeshType::New();

  typedef itk::MeshFunctionTestHelper<MeshType, double> InterpolatorType;

  InterpolatorType::Pointer interpolator = InterpolatorType::New();

  interpolator->SetInputMesh( mesh );

  interpolator->RunTest();

  std::cout << "Test passed" << std::endl;

  return EXIT_SUCCESS;
}
