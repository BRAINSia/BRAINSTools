/*=========================================================================
 *
 *  Copyright SINAPSE: Scalable Informatics for Neuroscience, Processing and Software Engineering
 *            The University of Iowa
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         http://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *=========================================================================*/
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

int
main( int argc, char * argv[] )
{
  PARSE_ARGS;
  BRAINSRegisterAlternateIO();

  if ( inputSurfaceFile == "" )
  {
    std::cerr << "No input file specified" << std::endl;
    return 1;
  }
  if ( outputSurfaceFile == "" )
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

  using MeshType = itk::QuadEdgeMesh< float, 3 >;

  // read original mesh
  using ReaderType = itk::QuadEdgeMeshVTKPolyDataReader< MeshType >;

  ReaderType::Pointer inputReader = ReaderType::New();
  inputReader->SetFileName( inputSurfaceFile.c_str() );
  inputReader->Update();

  // writer type
  using WriterType = itk::QuadEdgeMeshScalarDataVTKPolyDataWriter< MeshType >;

  WriterType::Pointer writer = WriterType::New();

  // define reample type
  using TransformType = itk::IdentityTransform< double >;

  TransformType::Pointer transform = TransformType::New();

  using LinearInterpolatorType = itk::LinearInterpolateMeshFunction< MeshType >;
  using NearestInterpolatorType = itk::NearestNeighborInterpolateMeshFunction< MeshType >;

  LinearInterpolatorType::Pointer interpolator_l = LinearInterpolatorType::New();
  interpolator_l->SetUseNearestNeighborInterpolationAsBackup( false );

  NearestInterpolatorType::Pointer interpolator_n = NearestInterpolatorType::New();

  using ResamplingFilterType = itk::ResampleQuadEdgeMeshFilter< MeshType, MeshType >;

  ResamplingFilterType::Pointer resampler = ResamplingFilterType::New();

  resampler->SetTransform( transform );

  // set the interpolation type
  if ( interpolateType == "Nearest" )
  {
    resampler->SetInterpolator( interpolator_n );
  }
  else if ( interpolateType == "Linear" )
  {
    resampler->SetInterpolator( interpolator_l );
  }

  // sphere source type
  using SphereMeshSourceType = itk::IcosahedralRegularSphereMeshSource< MeshType >;

  SphereMeshSourceType::Pointer sphereMeshSource = SphereMeshSourceType::New();

  using PointType = SphereMeshSourceType::PointType;
  using VectorType = PointType::VectorType;

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
  resampler->SetReferenceMesh( sphereMeshSource->GetOutput() );
  resampler->SetInput( inputReader->GetOutput() );
  resampler->Update();

  writer->SetInput( resampler->GetOutput() );
  writer->SetFileName( outputSurfaceFile.c_str() );
  writer->Update();

  return 0;
}
