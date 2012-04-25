/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkMeanSquaresMeshToMeshMetricTest1.cxx,v $
  Language:  C++
  Date:      $Date: 2007-09-06 17:44:24 $
  Version:   $Revision: 1.3 $

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
#include "itkRegularSphereMeshSource.h"
#include "itkQuadEdgeMeshScalarDataVTKPolyDataWriter.h"

int main( int argc, char * argv [] )
{
  if( argc < 4 )
    {
    std::cerr << "Missing arguments" << std::endl;
    std::cerr << "Usage: " << std::endl;
    std::cerr << argv[0] << std::endl;
    std::cerr << "outputSphericalMesh.vtk ";
    std::cerr << " radius (mm) ";
    std::cerr << "resolution(4=IC4,5=IC5...) " << std::endl;
    return EXIT_FAILURE;
    }

  typedef itk::QuadEdgeMesh<float, 3> MeshType;

  typedef itk::RegularSphereMeshSource<MeshType> SphereMeshSourceType;

  SphereMeshSourceType::Pointer sphereMeshSource = SphereMeshSourceType::New();

  typedef SphereMeshSourceType::PointType PointType;
  typedef PointType::VectorType           VectorType;

  PointType center;
  center.Fill( 0.0 );

  VectorType scaleVector;
  scaleVector.Fill( atof( argv[2] ) );

  const unsigned int resolution = atoi( argv[3] );

  sphereMeshSource->SetCenter( center );
  sphereMeshSource->SetScale( scaleVector );
  sphereMeshSource->SetResolution( resolution );

  try
    {
    sphereMeshSource->Update();
    }
  catch( itk::ExceptionObject & excp )
    {
    std::cerr << "Error during writer Update() " << std::endl;
    std::cerr << excp << std::endl;
    return EXIT_FAILURE;
    }

  typedef itk::QuadEdgeMeshScalarDataVTKPolyDataWriter<MeshType> WriterType;

  WriterType::Pointer writer = WriterType::New();

  writer->SetInput( sphereMeshSource->GetOutput()  );
  writer->SetFileName( argv[1] );

  try
    {
    writer->Update();
    }
  catch( itk::ExceptionObject & excp )
    {
    std::cerr << "Error during writer Update() " << std::endl;
    std::cerr << excp << std::endl;
    return EXIT_FAILURE;
    }

  return EXIT_SUCCESS;
}
