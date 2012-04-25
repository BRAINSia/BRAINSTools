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

#include "itkQuadEdgeMeshTraits.h"

#include "itkQuadEdgeMeshVTKPolyDataReader.h"
#include "itkQuadEdgeMeshScalarDataVTKPolyDataWriter.h"

#include "itkLinearInterpolateMeshFunction.h"
#include "itkWarpQuadEdgeMeshFilter.h"

int main( int argc, char * argv [] )
{
  if( argc < 5 )
    {
    std::cerr << "Missing arguments" << std::endl;
    std::cerr << "Usage: " << std::endl;
    std::cerr << argv[0] << std::endl;
    std::cerr << "inputFixedMesh inputReferenceMesh ";
    std::cerr << "inputDeformationField outputWarpedMesh ";
    std::cerr << std::endl;
    return EXIT_FAILURE;
    }

  typedef float MeshPixelType;
  const unsigned int Dimension = 3;

  typedef itk::QuadEdgeMesh<MeshPixelType, Dimension> InputMeshType;
  typedef itk::QuadEdgeMesh<MeshPixelType, Dimension> ReferenceMeshType;

  typedef InputMeshType::PointType PointType;
  typedef PointType::VectorType    VectorType;

  typedef itk::QuadEdgeMeshTraits<VectorType, Dimension, bool, bool>     VectorPointSetTraits;
  typedef itk::QuadEdgeMesh<VectorType, Dimension, VectorPointSetTraits> MeshWithVectorsType;

  typedef itk::QuadEdgeMeshVTKPolyDataReader<MeshWithVectorsType> DeformationFieldReaderType;
  typedef itk::QuadEdgeMeshVTKPolyDataReader<InputMeshType>       InputReaderType;
  typedef itk::QuadEdgeMeshVTKPolyDataReader<ReferenceMeshType>   ReferenceReaderType;

  InputReaderType::Pointer            inputMeshReader = InputReaderType::New();
  ReferenceReaderType::Pointer        referenceMeshReader = ReferenceReaderType::New();
  DeformationFieldReaderType::Pointer deformationFieldReader
    = DeformationFieldReaderType::New();

  inputMeshReader->SetFileName( argv[1] );
  referenceMeshReader->SetFileName( argv[2] );
  deformationFieldReader->SetFileName( argv[3] );

  typedef itk::LinearInterpolateMeshFunction<ReferenceMeshType> InterpolatorType;

  InterpolatorType::Pointer interpolator = InterpolatorType::New();

  typedef itk::WarpQuadEdgeMeshFilter<
      InputMeshType, ReferenceMeshType, MeshWithVectorsType> WarpFilterType;

  WarpFilterType::Pointer warpFilter = WarpFilterType::New();

  warpFilter->SetInterpolator(interpolator);

  warpFilter->SetInputMesh(inputMeshReader->GetOutput() );

  warpFilter->SetReferenceMesh(referenceMeshReader->GetOutput() );

  warpFilter->SetDeformationField(deformationFieldReader->GetOutput() );

  warpFilter->Update();

  // write the result
  typedef itk::QuadEdgeMeshScalarDataVTKPolyDataWriter<InputMeshType> WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetInput( warpFilter->GetOutput() );
  writer->SetFileName(argv[4]);
  writer->Update();

  return EXIT_SUCCESS;
}
