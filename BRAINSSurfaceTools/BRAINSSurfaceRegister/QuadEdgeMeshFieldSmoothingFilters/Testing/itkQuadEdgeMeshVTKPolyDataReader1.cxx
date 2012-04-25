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

#include "itkQuadEdgeMesh.h"
#include "itkMeshWriterHelper1.h"
#include "itkQuadEdgeMeshVTKPolyDataReader.h"
#include "itkTestingMacros.h"
#include "itksys/SystemTools.hxx"

int main( int argc, char *argv[] )
{
  if( argc < 3 )
    {
    std::cerr << "Missing Parameters " << std::endl;
    std::cerr << "Usage: " << argv[0];
    std::cerr << " inputMeshFile  outputMeshfile " << std::endl;
    return EXIT_FAILURE;
    }

  const unsigned int Dimension = 3;
  typedef float MeshPixelType;

  typedef itk::QuadEdgeMesh<MeshPixelType, Dimension> MeshType;

  typedef itk::QuadEdgeMeshVTKPolyDataReader<MeshType> ReaderType;

  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[1] );

  try
    {
    reader->Update();
    }
  catch( itk::ExceptionObject & exp )
    {
    std::cerr << "Exception thrown while reading the input file " << std::endl;
    std::cerr << exp << std::endl;
    return EXIT_FAILURE;
    }

  itk::MeshWriterHelper1<MeshType>::WriteMeshToFile( reader->GetOutput(), argv[2] );

  if( itksys::SystemTools::FilesDiffer( argv[1], argv[2] ) )
    {
    std::cerr << "Output file differs from input file" << std::endl;
    std::cerr << "This may indicate an error in the Reader or the Writer" << std::endl;
    return EXIT_FAILURE;
    }

  return EXIT_SUCCESS;
}
