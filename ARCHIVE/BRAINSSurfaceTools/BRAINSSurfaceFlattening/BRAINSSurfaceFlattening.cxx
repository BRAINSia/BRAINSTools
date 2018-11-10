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
 Date:      $Date: 2011/09/09 14:53:40 $
 Version:   $Revision: 1.0 $

   Copyright (c) University of Iowa Department of Radiology. All rights reserved.
   See GTRACT-Copyright.txt or http://mri.radiology.uiowa.edu/copyright/GTRACT-Copyright.txt
   for details.

      This software is distributed WITHOUT ANY WARRANTY; without even
      the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
      PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#include "itkQuadEdgeMesh.h"
#include "itkTransformMeshFilter.h"
#include "itkVersorTransform.h"

#include "itkVTKPolyDataReader.h"
#include "itkVTKPolyDataWriter.h"

#include "VNLIterativeSparseSolverTraits.h"
#include "itkQuadEdgeMeshParamMatrixCoefficients.h"

#include "itkQuadEdgeMeshToSphereFilter.h"
#include "itkVector.h"

#include "BRAINSSurfaceFlatteningCLP.h"
#include <BRAINSCommonLib.h>

int main( int argc, char* * argv )
{
  // if( argc != 7 )
  //   {
  //   std::cout <<" It requires 6 arguments " <<std::endl;
  //   std::cout <<" InputFileName OutputFileName " <<std::endl;
  //   std::cout <<" radius seedAxis(0|1|2) " <<std::endl;
  // std::cout <<" maxAxis(0|1|2) alignAxis(0|1|2) " <<std::endl;
  //   return EXIT_FAILURE;
  //   }

  PARSE_ARGS;
  BRAINSRegisterAlternateIO();
  // check input file name
  if( inputSurfaceFile == "" )
    {
    std::cerr << "No input surface file specified" << std::endl;
    return EXIT_FAILURE;
    }
  // check output file name
  if( outputSurfaceFile == "" )
    {
    std::cerr << "No output surface file specified" << std::endl;
    return EXIT_FAILURE;
    }
  std::cout << "Surface Flattening Parameters" << std::endl;
  std::cout << "------------------------------------------------------" << std::endl;
  std::cout << "\tInput Surface FileName: " << inputSurfaceFile << std::endl;
  std::cout << "\tOutput Surface FileName: " << outputSurfaceFile << std::endl;
  std::cout << "\tThe radius of the output sphere: " << sphereRadius << std::endl;
  if( seed )
    {
    std::cout << "\tSpecify the seeds along " << seedAxis << std::endl;
    }
  else
    {
    std::cout << "\tSepcify the seeds using max and min cell Ids." << std::endl;
    }
  if( rotate )
    {
    std::cout << "\tRotate the sphere so that polar point along " << maxAxis;
    std::cout << " remains at the top." << std::endl;
    }
  else
    {
    std::cout << "\tDo not rotate the sphere." << std::endl;
    }
  std::cout << "------------------------------------------------------" << std::endl;

  // ** TYPEDEF **
  using Coord = double;

  using MeshType = itk::QuadEdgeMesh<Coord, 3>;
  using MeshPointer = MeshType::Pointer;
  using CellIdentifier = MeshType::CellIdentifier;
  using PointIdentifier = MeshType::PointIdentifier;
  using QEPrimal = MeshType::QEPrimal;

  using PointsContainerConstIterator = MeshType::PointsContainerConstIterator;

  using ReaderType = itk::VTKPolyDataReader<MeshType>;
  using WriterType = itk::VTKPolyDataWriter<MeshType>;

  using SolverTraits = VNLIterativeSparseSolverTraits<Coord>;

  using FilterType = itk::QuadEdgeMeshToSphereFilter<
      MeshType, MeshType, SolverTraits>;

  using TransformType = itk::VersorTransform<Coord>;
  using VersorType = TransformType::VersorType;
  using TransformMeshType = itk::TransformMeshFilter<MeshType, MeshType, TransformType>;

  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( inputSurfaceFile.c_str() );
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

  MeshPointer mesh = reader->GetOutput();

  // choose the seeds along the specified seedAxis
  unsigned short seed_ax = 0;
  if( seedAxis == "y" )
    {
    seed_ax = 1;
    }
  else if( seedAxis == "z" )
    {
    seed_ax = 2;
    }

  // fix the polar point along the specified maxAxis
  unsigned short max_ax = 0;
  if( maxAxis == "y" )
    {
    max_ax = 1;
    }
  else if( maxAxis == "z" )
    {
    max_ax = 2;
    }

  // walk through all of the points of the input surface
  // get the one with max and min coordinates along the seedAxis
  float seedMin = 999999.0;
  float seedMax = -999999.0;

  // get the one with max coordinates along the maxAxis
  float polarMax = -999999.0;

  MeshType::PointsContainerConstPointer points = mesh->GetPoints();

  PointsContainerConstIterator p_It = points->begin();

  CellIdentifier seedCell0, seedCell1;

  PointIdentifier polarPId = 0;

  while( p_It != points->end() )
    {
    MeshType::PointType p = p_It.Value();

    // look for seedCell0 along seedAxis-
    if( p[seed_ax] < seedMin )
      {
      if( p.GetEdge() != (QEPrimal *)nullptr )
        {
        QEPrimal* e = p.GetEdge();
        if( e->GetLeft() != mesh->m_NoFace )
          {
          seedCell0 = e->GetLeft();
          seedMin = p[seed_ax];
          }
        }
      }

    // look for seedCell1 along seedAxis+
    if( p[seed_ax] > seedMax )
      {
      if( p.GetEdge() != (QEPrimal *)nullptr )
        {
        QEPrimal* e = p.GetEdge();
        if( e->GetLeft() != mesh->m_NoFace )
          {
          seedCell1 = e->GetLeft();
          seedMax = p[seed_ax];
          }
        }
      }

    // look for the polar point along maxAxis+
    if( p[max_ax] > polarMax )
      {
      polarMax = p[max_ax];
      polarPId = p_It.Index();
      }

    p_It++;
    }

  std::vector<CellIdentifier> seedList;

  if( seed ) // user specify seeds
    {
    seedList.push_back(seedCell0);
    seedList.push_back(seedCell1);
    }
  else
    {
    using CellsContainerConstPointer = MeshType::CellsContainerConstPointer;
    CellsContainerConstPointer cells = mesh->GetCells();

    using CellsContainerConstIterator = MeshType::CellsContainerConstIterator;
    CellsContainerConstIterator c_It = cells->begin();

    seedList.push_back(c_It.Index() );
    c_It = cells->end();
    --c_It;

    seedList.push_back(c_It.Index() );
    }

  itk::OnesMatrixCoefficients<MeshType> coeff0;

  FilterType::Pointer filter = FilterType::New();
  filter->SetInput( mesh );
  filter->SetCoefficientsMethod( &coeff0 );
  filter->SetRadius(sphereRadius);
  filter->SetSeedFaces(seedList);
  filter->Update();

  MeshPointer sphere0 = filter->GetOutput();

  // rotate sphere0 so that
  // the polar point remains at the polar along the maxAxis
  // ***QuadEdgeMeshToSphereFilter generates the sphere
  // ***with seedCells along the z axis.
  if( rotate )
    {
    TransformMeshType::Pointer transformFilter = TransformMeshType::New();
    TransformType::Pointer     transform = TransformType::New();
    transform->SetIdentity();

    VersorType::VectorType axis;
    axis.Fill(0.0);
    unsigned short ax = 2; // always z axis
    axis[ax] = 1.0;

    // caculate the angle
    MeshType::PointType polarPOnSphere = sphere0->GetPoint(polarPId);
    using NewVectorType = itk::Vector<Coord, 2>;
    NewVectorType polarPOnSphere2D;
    polarPOnSphere2D[0] = polarPOnSphere[0];
    polarPOnSphere2D[1] = polarPOnSphere[1];
    if( polarPOnSphere2D[1] > 0.0 )
      {
      axis[ax] = -1.0;
      }

    NewVectorType alignVector2D;
    alignVector2D.Fill(0.0);
    // align_ax should be either x or y
    unsigned short align_ax = 0; // x
    alignVector2D[align_ax] = 1.0;

    VersorType::ValueType angle;
    angle = acos( (polarPOnSphere2D * alignVector2D) / (polarPOnSphere2D.GetNorm() * alignVector2D.GetNorm() ) );

    VersorType rotation;
    rotation.Set( axis, angle );
    transform->SetRotation( rotation );

    transformFilter->SetInput(sphere0);
    transformFilter->SetTransform(transform);
    transformFilter->Update();

    WriterType::Pointer writer = WriterType::New();
    writer->SetInput( transformFilter->GetOutput() );
    writer->SetFileName( outputSurfaceFile.c_str() );
    writer->Update();
    }
  else
    {
    WriterType::Pointer writer = WriterType::New();
    writer->SetInput( sphere0 );
    writer->SetFileName( outputSurfaceFile.c_str() );
    writer->Update();
    }

  return EXIT_SUCCESS;
}
