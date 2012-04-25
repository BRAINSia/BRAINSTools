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

#include "itkMeshGeneratorHelper5.h"
#include "itkMeshWriterHelper1.h"

int main( int argc, char *argv[] )
{
  if( argc < 5 )
    {
    std::cerr << "Missing Parameters " << std::endl;
    std::cerr << "Usage: " << argv[0];
    std::cerr << " outputMeshFile resolution scale factor ";
    return EXIT_FAILURE;
    }

  const unsigned int Dimension = 3;
  typedef float MeshPixelType;

  typedef itk::QuadEdgeMesh<MeshPixelType, Dimension> MeshType;

  typedef itk::MeshGeneratorHelper5<MeshType> GeneratorType;

  MeshType::Pointer mesh;

  const unsigned int resolution = atoi( argv[2] );
  const double       scale = atof( argv[3] );
  const double       factor = atof( argv[4] );

  GeneratorType::GenerateMesh( mesh, resolution, scale, factor );

  itk::MeshWriterHelper1<MeshType>::WriteMeshToFile( mesh, argv[1] );

  return EXIT_SUCCESS;
}
