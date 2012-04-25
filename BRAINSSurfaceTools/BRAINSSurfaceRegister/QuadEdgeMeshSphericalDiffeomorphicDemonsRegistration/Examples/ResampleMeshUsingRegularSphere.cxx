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

#include "itkLinearInterpolateMeshFunction.h"
#include "itkIdentityTransform.h"
#include "itkQuadEdgeMesh.h"

#include "itkCommand.h"
#include "itkVTKPolyDataReader.h"

#include "itkResampleQuadEdgeMeshFilter.h"
#include "itkQuadEdgeMeshScalarDataVTKPolyDataWriter.h"

int main( int argc, char * argv [] )
{
  if( argc < 3 )
    {
    std::cerr << "Missing arguments" << std::endl;
    std::cerr << "Usage: " << std::endl;
    std::cerr << argv[0] << std::endl;
    std::cerr << "inputFixedMesh inputMovingMesh ";
    std::cerr << "outputResampledMesh " << std::endl;
    return EXIT_FAILURE;
    }

  typedef itk::QuadEdgeMesh<float, 3> FixedMeshType;
  typedef itk::QuadEdgeMesh<float, 3> MovingMeshType;

  typedef itk::VTKPolyDataReader<FixedMeshType>  FixedReaderType;
  typedef itk::VTKPolyDataReader<MovingMeshType> MovingReaderType;

  FixedReaderType::Pointer fixedMeshReader = FixedReaderType::New();
  fixedMeshReader->SetFileName( argv[1] );

  MovingReaderType::Pointer movingMeshReader = MovingReaderType::New();
  movingMeshReader->SetFileName( argv[2] );

  try
    {
    fixedMeshReader->Update();
    movingMeshReader->Update();
    }
  catch( itk::ExceptionObject & exp )
    {
    std::cerr << exp << std::endl;
    return EXIT_FAILURE;
    }

  FixedMeshType::ConstPointer  meshFixed  = fixedMeshReader->GetOutput();
  MovingMeshType::ConstPointer meshMoving = movingMeshReader->GetOutput();

  typedef itk::IdentityTransform<double> TransformType;

  TransformType::Pointer transform = TransformType::New();

  typedef itk::LinearInterpolateMeshFunction<MovingMeshType> InterpolatorType;

  InterpolatorType::Pointer interpolator = InterpolatorType::New();
  interpolator->SetUseNearestNeighborInterpolationAsBackup(false);
  std::cout << interpolator->GetUseNearestNeighborInterpolationAsBackup() << std::endl;

  typedef itk::ResampleQuadEdgeMeshFilter<FixedMeshType, MovingMeshType> ResamplingFilterType;

  ResamplingFilterType::Pointer resampler = ResamplingFilterType::New();

  resampler->SetReferenceMesh( meshFixed );
  resampler->SetInput( meshMoving );

  resampler->SetTransform( transform );
  resampler->SetInterpolator( interpolator );

  resampler->Update();

  typedef itk::QuadEdgeMeshScalarDataVTKPolyDataWriter<FixedMeshType> WriterType;

  WriterType::Pointer writer = WriterType::New();

  writer->SetInput( resampler->GetOutput()  );
  writer->SetFileName( argv[3] );

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
