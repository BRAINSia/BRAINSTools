/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkKdTreeGeneratorTest.cxx,v $
  Language:  C++
  Date:      $Date: 2005-07-26 15:55:12 $
  Version:   $Revision: 1.4 $

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
#include "itkVTKPolyDataReader.h"
#include "itkTimeProbesCollectorBase.h"
#include "itkPointLocator2.h"

int main(int argc, char *argv[])
{
  if( argc < 2 )
    {
    std::cerr << "Usage: " << argv[0] << " vtkInputFilename " << std::endl;
    return EXIT_FAILURE;
    }

  const unsigned int Dimension = 3;

  typedef itk::QuadEdgeMesh<double, Dimension> MeshType;

  typedef itk::VTKPolyDataReader<MeshType> ReaderType;

  itk::TimeProbesCollectorBase chronometer;

  // Read mesh file

  std::cout << "Read " << argv[1] << std::endl;

  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[1] );

  chronometer.Start("Reading");

  try
    {
    reader->Update();
    }
  catch( itk::ExceptionObject & excp )
    {
    std::cerr << excp << std::endl;
    return EXIT_FAILURE;
    }

  chronometer.Stop("Reading");

  MeshType::ConstPointer mesh = reader->GetOutput();

  // Instantiate the point locator
  typedef itk::PointLocator2<MeshType> PointLocatorType;

  PointLocatorType::Pointer locator = PointLocatorType::New();

  locator->SetPointSet( mesh );

  MeshType::ConstPointer mesh2 = locator->GetPointSet();

  if( mesh.GetPointer() != mesh2.GetPointer() )
    {
    std::cerr << "Error in SetPointSet()/GetPointSet() " << std::endl;
    return EXIT_FAILURE;
    }

  chronometer.Start("Generate");

  try
    {
    locator->Initialize();
    }
  catch( itk::ExceptionObject & excp )
    {
    std::cerr << excp << std::endl;
    return EXIT_FAILURE;
    }

  chronometer.Stop("Generate");

  std::cout << locator->GetNameOfClass() << std::endl;
  locator->Print( std::cout );

  unsigned int numberOfPoints = mesh->GetNumberOfPoints();

  unsigned int numberOfNeighbors = 1;

  PointLocatorType::InstanceIdentifierVectorType neighbors(numberOfNeighbors);

  MeshType::PointType queryPoint;

  bool testPassed = true;

  chronometer.Start("SearchClosest");
  for( unsigned int i = 0; i < numberOfPoints; i++ )
    {
    mesh->GetPoint(i, &queryPoint);

    locator->Search( queryPoint, numberOfNeighbors, neighbors );

    if( neighbors[0] != i )
      {
      std::cerr << "difference at " << i << " " << neighbors[0] << std::endl;
      testPassed = false;
      }
    }

  chronometer.Stop("SearchClosest");

  chronometer.Report( std::cout );

  if( !testPassed )
    {
    std::cerr << "Test Failed." << std::endl;
    return EXIT_FAILURE;
    }

  std::cout << "Test passed." << std::endl;

  return EXIT_SUCCESS;
}
