/*=========================================================================

 Program:   BRAINS (Brain Research: Analysis of Images, Networks, and Systems)
 Module:    $RCSfile: $
 Language:  C++
 Date:      $Date: 2011/07/09 14:53:40 $
 Version:   $Revision: 1.0 $

 Copyright (c) University of Iowa Department of Radiology. All rights reserved.
 See GTRACT-Copyright.txt or http://mri.radiology.uiowa.edu/copyright/GTRACT-Copyright.txt
 for details.

 This software is distributed WITHOUT ANY WARRANTY; without even
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE.  See the above copyright notices for more information.

 =========================================================================*/

#include "itkLinearInterpolateMeshFunction.h"
#include "itkNearestNeighborInterpolateMeshFunction.h"
#include "itkIdentityTransform.h"
#include "itkQuadEdgeMesh.h"
#include "itkIcosahedralRegularSphereMeshSource.h"

#include "itkQuadEdgeMeshScalarPixelValuesSmoothingFilter.h"
#include "itkResampleQuadEdgeMeshFilter.h"
#include "itkRescaleScalarsQuadEdgeMeshFilter.h"

#include "itkQuadEdgeMeshVTKPolyDataReader.h"
#include "itkQuadEdgeMeshScalarDataVTKPolyDataWriter.h"

#include "IcosahedronResamplerCLP.h"
#include <BRAINSCommonLib.h>

int main( int argc, char * argv [] )
{
  PARSE_ARGS;
  BRAINSRegisterAlternateIO();

  if( inputSurfaceFile == "" )
    {
    std::cerr << "No input file specified" << std::endl;
    return 1;
    }
  if( outputSurfaceFile == "" )
    {
    std::cerr << "No output file specified" << std::endl;
    return 1;
    }

  std::cout << "---------------------------------------------------" << std::endl;
  std::cout << "Input Surface: " << inputSurfaceFile << std::endl;
  std::cout << "Output Surface: " << outputSurfaceFile << std::endl;
  std::cout << "Resample the input surface with an icosahedron sphere ";
  std::cout << "at the resolution level of " << resolution;
  std::cout << " with the radius of " << radius << std::endl;
  std::cout << "Interpolation Type: " << interpolateType << std::endl;
  std::cout << "---------------------------------------------------" << std::endl;

  typedef itk::QuadEdgeMesh<float, 3> MeshType;

  // read original mesh
  typedef itk::QuadEdgeMeshVTKPolyDataReader<MeshType> ReaderType;

  ReaderType::Pointer inputReader = ReaderType::New();
  inputReader->SetFileName( inputSurfaceFile.c_str() );
  inputReader->Update();

  // writer type
  typedef itk::QuadEdgeMeshScalarDataVTKPolyDataWriter<MeshType> WriterType;

  WriterType::Pointer writer = WriterType::New();

  // define reample type
  typedef itk::IdentityTransform<double> TransformType;

  TransformType::Pointer transform = TransformType::New();

  typedef itk::LinearInterpolateMeshFunction<MeshType>          LinearInterpolatorType;
  typedef itk::NearestNeighborInterpolateMeshFunction<MeshType> NearestInterpolatorType;

  LinearInterpolatorType::Pointer interpolator_l = LinearInterpolatorType::New();
  interpolator_l->SetUseNearestNeighborInterpolationAsBackup(false);

  NearestInterpolatorType::Pointer interpolator_n = NearestInterpolatorType::New();

  typedef itk::ResampleQuadEdgeMeshFilter<MeshType, MeshType> ResamplingFilterType;

  ResamplingFilterType::Pointer resampler = ResamplingFilterType::New();

  resampler->SetTransform( transform );

  // set the interpolation type
  if( interpolateType == "Nearest" )
    {
    resampler->SetInterpolator( interpolator_n );
    }
  else if( interpolateType == "Linear" )
    {
    resampler->SetInterpolator( interpolator_l );
    }

  // sphere source type
  typedef itk::IcosahedralRegularSphereMeshSource<MeshType> SphereMeshSourceType;

  SphereMeshSourceType::Pointer sphereMeshSource = SphereMeshSourceType::New();

  typedef SphereMeshSourceType::PointType PointType;
  typedef PointType::VectorType           VectorType;

  PointType center;
  center.Fill( 0.0 );

  VectorType scaleVector;
  scaleVector.Fill( radius );

  // generate the sphere source
  sphereMeshSource->SetCenter( center );
  sphereMeshSource->SetScale( scaleVector );
  sphereMeshSource->SetResolution( resolution );
  sphereMeshSource->Update();

  // resample
  resampler->SetReferenceMesh(sphereMeshSource->GetOutput() );
  resampler->SetInput( inputReader->GetOutput() );
  resampler->Update();

  writer->SetInput(resampler->GetOutput() );
  writer->SetFileName(outputSurfaceFile.c_str() );
  writer->Update();

  return 0;
}
