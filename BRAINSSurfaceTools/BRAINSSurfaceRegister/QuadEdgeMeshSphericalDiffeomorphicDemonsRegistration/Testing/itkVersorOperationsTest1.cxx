/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: ImageRegistration8.cxx,v $
  Language:  C++
  Date:      $Date: 2008-03-10 19:46:31 $
  Version:   $Revision: 1.37 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif

#include "itkVersor.h"

int main( int argc, char *argv[] )
{
  typedef double ComponentType;

  typedef itk::Versor<ComponentType> VersorType;
  typedef VersorType::VectorType     AxisType;
  typedef VersorType::ValueType      AngleType;

  AxisType axis1;

  axis1[0] = 1.0;
  axis1[1] = 0.0;
  axis1[2] = 0.0;

  AngleType angle1 = vnl_math::pi / 2.0;

  VersorType versor1;
  versor1.Set( axis1, angle1 );

  std::cout << "versor1 = " << versor1 << std::endl;

  AxisType axis2;

  axis2[0] = 0.0;
  axis2[1] = 0.0;
  axis2[2] = 1.0;

  AngleType angle2 = vnl_math::pi / 2.0;

  VersorType versor2;
  versor2.Set( axis2, angle2 );

  std::cout << "versor2 = " << versor2 << std::endl;

  VersorType versor12 = versor2 * versor1;

  std::cout << "versor12 = " << versor12 << std::endl;

  return EXIT_SUCCESS;
}
