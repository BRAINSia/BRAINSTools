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
#include "itkVector.h"
#include "itkListSample.h"
#include "itkKdTree.h"
#include "itkKdTreeGenerator.h"
#include "itkPointSetToListSampleAdaptor.h"
#include "itkTimeProbesCollectorBase.h"

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

  MeshType::Pointer mesh = reader->GetOutput();

  // Adapt the Mesh to make it look like a ListSample
  typedef itk::Statistics::PointSetToListSampleAdaptor<MeshType> SampleType;

  SampleType::Pointer sample = SampleType::New();

  sample->SetMeasurementVectorSize( Dimension );

  sample->SetPointSet( mesh );

  // This is the same PointType of the MeshType;
  typedef SampleType::MeasurementVectorType MeasurementVectorType;

  // Instantiate the KdTreeGenerator
  typedef itk::Statistics::KdTreeGenerator<SampleType> TreeGeneratorType;
  TreeGeneratorType::Pointer treeGenerator = TreeGeneratorType::New();

  treeGenerator->SetSample( sample );
  treeGenerator->SetBucketSize( 16 );

  chronometer.Start("Generate");

  try
    {
    treeGenerator->Update();
    }
  catch( itk::ExceptionObject & excp )
    {
    std::cerr << excp << std::endl;
    return EXIT_FAILURE;
    }

  chronometer.Stop("Generate");

  typedef TreeGeneratorType::KdTreeType TreeType;
  typedef TreeType::NearestNeighbors    NeighborsType;
  typedef TreeType::KdTreeNodeType      NodeType;

  TreeType::Pointer tree = treeGenerator->GetOutput();

  typedef MeasurementVectorType PointType;

  unsigned int numberOfPoints = mesh->GetNumberOfPoints();

  TreeType::InstanceIdentifierVectorType neighbors;
  unsigned int                           numberOfNeighbors = 1;

  PointType queryPoint;

  bool testPassed = true;

  chronometer.Start("SearchClosest");
  for( unsigned int i = 0; i < numberOfPoints; i++ )
    {
    mesh->GetPoint(i, &queryPoint);
    tree->Search( queryPoint, numberOfNeighbors, neighbors );

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
