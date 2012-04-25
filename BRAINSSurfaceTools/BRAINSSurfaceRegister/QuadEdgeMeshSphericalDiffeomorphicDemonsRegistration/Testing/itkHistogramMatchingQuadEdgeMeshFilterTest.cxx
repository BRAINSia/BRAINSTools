/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkHistogramMatchingQuadEdgeMeshFilterTest.cxx,v $
  Language:  C++
  Date:      $Date: 2010-05-17 19:46:31 $
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
#include "itkQuadEdgeMeshVTKPolyDataReader.h"
#include "itkQuadEdgeMeshScalarDataVTKPolyDataWriter.h"

#include "itkQuadEdgeMesh.h"
#include "itkHistogramMatchingQuadEdgeMeshFilter.h"

int main( int argc, char * argv [] )
{
  if( argc < 4 )
    {
    std::cerr << "Missing arguments" << std::endl;
    std::cerr << "Usage: " << std::endl;
    std::cerr << argv[0] << std::endl;
    std::cerr << "inputPolyData referencePolyData ";
    std::cerr << "outputPolyData";
    std::cerr << "NumberOfHistogramLevels NumberOfMatchPoints";
    std::cerr << std::endl;
    return EXIT_FAILURE;
    }

  typedef float MeshPixelType;
  const unsigned int Dimension = 3;

  typedef itk::QuadEdgeMesh<MeshPixelType, Dimension> MeshType;

  typedef itk::QuadEdgeMeshVTKPolyDataReader<MeshType> SourceReaderType;
  typedef itk::QuadEdgeMeshVTKPolyDataReader<MeshType> ReferenceReaderType;

  SourceReaderType::Pointer srcReader = SourceReaderType::New();
  srcReader->SetFileName( argv[1] );

  ReferenceReaderType::Pointer refReader = ReferenceReaderType::New();
  refReader->SetFileName( argv[2] );

  srcReader->Update();
  refReader->Update();

  typedef itk::HistogramMatchingQuadEdgeMeshFilter<MeshType, MeshType> FilterType;

  FilterType::Pointer filter  = FilterType::New();

  filter->SetSourceMesh(srcReader->GetOutput() );
  filter->SetReferenceMesh(refReader->GetOutput() );
  filter->SetNumberOfHistogramLevels( atoi(argv[4]) );
  filter->SetNumberOfMatchPoints( atoi(argv[5]) );
  filter->Update();

  // Walk the output and compare with reference
  MeshType::Pointer output = filter->GetOutput();
  MeshType::Pointer ref    = refReader->GetOutput();

  typedef MeshType::PointDataContainerConstPointer PointDataContainerConstPointer;
  PointDataContainerConstPointer output_pointData = output->GetPointData();
  PointDataContainerConstPointer ref_pointData = ref->GetPointData();

  typedef MeshType::PointDataContainer::ConstIterator PointDataIterator;
  PointDataIterator out_Itr = output_pointData->Begin();
  PointDataIterator out_End = output_pointData->End();

  PointDataIterator ref_Itr = ref_pointData->Begin();
  PointDataIterator ref_End = ref_pointData->End();

  bool passed = true;
  while( out_Itr != out_End )
    {
    MeshType::PixelType diff = ref_Itr.Value() - out_Itr.Value();
    if( vnl_math_abs( diff ) > 0.1 )
      {
      passed = false;
      std::cout << "Test failed at: " << out_Itr.Index() << " ";
      std::cout << "Output value: " << out_Itr.Value() << " ";
      std::cout << "Ref value: " << ref_Itr.Value() << std::endl;
      }
    ++out_Itr;
    ++ref_Itr;
    }

  if( !passed )
    {
    std::cout << "Test failed." << std::endl;
    return EXIT_FAILURE;
    }

  std::cout << "Test passed." << std::endl;

  typedef itk::QuadEdgeMeshScalarDataVTKPolyDataWriter<MeshType> WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetInput( filter->GetOutput() );
  writer->SetFileName(argv[3]);
  writer->Update();

  return EXIT_SUCCESS;
}
