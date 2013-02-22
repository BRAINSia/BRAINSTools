/*=========================================================================

 Program:   BRAINS (Brain Research: Analysis of Images, Networks, and Systems)
 Module:    $RCSfile: $
 Language:  C++
 Date:      $Date: 2011/07/01 14:53:40 $
 Version:   $Revision: 1.9 $

   Copyright (c) University of Iowa Department of Radiology. All rights reserved.
   See GTRACT-Copyright.txt or http://mri.radiology.uiowa.edu/copyright/GTRACT-Copyright.txt
   for details.

      This software is distributed WITHOUT ANY WARRANTY; without even
      the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
      PURPOSE.  See the above copyright notices for more information.

=========================================================================*/

#include "vtkITKArchetypeImageSeriesReader.h"
#include "vtkITKArchetypeImageSeriesScalarReader.h"

#include "vtkImageChangeInformation.h"
#include "vtkTransform.h"
#include <vtkImageData.h>

#include "vtkImageMarchingCubes.h"
#include "vtkPolyDataConnectivityFilter.h"
#include "vtkCleanPolyData.h"
#include "vtkPolyData.h"
#include "vtkPolyDataReader.h"
#include "vtkTransformPolyDataFilter.h"
#include "vtkExtractEdges.h"

#include "itkConvertVTKToQuadEdgeMeshFilter.h"
#include "itkQuadEdgeMesh.h"
#include "itkQuadEdgeMeshDecimationCriteria.h"
#include "itkSquaredEdgeLengthDecimationQuadEdgeMeshFilter.h"
#include "itkQuadEdgeMeshParamMatrixCoefficients.h"
#include "itkSmoothingQuadEdgeMeshFilter.h"

#include "itkVTKPolyDataWriter.h"

#include "BRAINSSurfaceGenerationCLP.h"

#include <vtksys/SystemTools.hxx>
#include <vtkSmartPointer.h>

int main( int argc, char * *argv )
{
  PARSE_ARGS;
  if( (inputImageFile == "") && (inputSurfaceFile == "") )
    {
    std::cerr << "No input file specified" << std::endl;
    return EXIT_FAILURE;
    }
  if( outputSurface == "" )
    {
    std::cerr << "No output surface file specified" << std::endl;
    return EXIT_FAILURE;
    }
  std::cout << "Surface Generation Parameters" << std::endl;
  std::cout << "------------------------------------------------------" << std::endl;
  if( inputImageFile != "" )
    {
    std::cout << "\tInput Image FileName: " << inputImageFile << std::endl;
    std::cout << "\tSurface Value: " << surfaceValue << std::endl;
    }

  if( inputSurfaceFile != "" )
    {
    std::cout << "\tInput Surface FileName: " << inputSurfaceFile << std::endl;
    }

  std::cout << "\tOutput Surface: " << outputSurface << std::endl;

  std::cout << "\tDecimate Flag: " << decimateSurface << std::endl;
  if( decimateSurface )
    {
    std::cout << "\tDecimate the Surface with Parameters: " << std::endl;
    std::cout << "\tNumber of Faces: " << numberOfElements << std::endl;
    }

  std::cout << "\tSmooth Flag: " << smoothSurface << std::endl;
  if( smoothSurface )
    {
    std::cout << "\tSmooth the Surface with Parameters: " << std::endl;
    std::cout << "\tNumber Of Iterations: " << numIterations << std::endl;
    std::cout << "\tRelaxation Factor: " << relaxationFactor << std::endl;
    }

  std::cout << "\tGenus Flag: " << genusNumber << std::endl;
  std::cout << "------------------------------------------------------" << std::endl;

  // Read InputFile

  vtkSmartPointer<vtkPolyData> surface = vtkSmartPointer<vtkPolyData>::New();

  if( inputImageFile != "" )
    {
    std::cout << "Generate a surface from the input image." << std::endl;
    // read imput image
    vtkSmartPointer<vtkITKArchetypeImageSeriesScalarReader> reader =
      vtkSmartPointer<vtkITKArchetypeImageSeriesScalarReader>::New();

    reader->SetArchetype(inputImageFile.c_str() );
    reader->SetOutputScalarTypeToDouble();
    reader->SetDesiredCoordinateOrientationToNative();
    reader->SetUseNativeOriginOn();
    reader->Update();

    vtkSmartPointer<vtkImageChangeInformation> ici =
      vtkSmartPointer<vtkImageChangeInformation>::New();
    ici->SetInput(reader->GetOutput() );
    ici->SetOutputSpacing( 1, 1, 1 );
    ici->SetOutputOrigin( 0, 0, 0 );
    ici->Update();

    vtkTransform * transformIJKtoRAS = vtkTransform::New();
    transformIJKtoRAS->SetMatrix(reader->GetRasToIjkMatrix() );
    transformIJKtoRAS->Inverse();

    vtkSmartPointer<vtkImageMarchingCubes> marchingcubes =
      vtkSmartPointer<vtkImageMarchingCubes>::New();
    marchingcubes->SetInput(ici->GetOutput() );
    marchingcubes->SetValue(0, surfaceValue);
    marchingcubes->ComputeScalarsOff();
    marchingcubes->ComputeNormalsOff();
    marchingcubes->ComputeGradientsOff();
    marchingcubes->Update();

    // largest connected region
    vtkSmartPointer<vtkPolyDataConnectivityFilter> largest =
      vtkSmartPointer<vtkPolyDataConnectivityFilter>::New();
    largest->SetInput(marchingcubes->GetOutput() );
    largest->SetExtractionModeToLargestRegion();
    largest->Update();

    vtkSmartPointer<vtkCleanPolyData> clean = vtkSmartPointer<vtkCleanPolyData>::New();
    clean->SetInput(largest->GetOutput() );
    clean->ConvertPolysToLinesOff();
    clean->ConvertLinesToPointsOff();
    clean->Update();

    vtkSmartPointer<vtkTransformPolyDataFilter> transformer =
      vtkSmartPointer<vtkTransformPolyDataFilter>::New();
    transformer->SetInput(clean->GetOutput() );
    transformer->SetTransform(transformIJKtoRAS);
    transformer->Update();

    surface = transformer->GetOutput();
    }

  if( inputSurfaceFile != "" )
    {
    std::cout << "Input surface: " << inputSurfaceFile.c_str() << std::endl;
    vtkSmartPointer<vtkPolyDataReader> surfaceReader = vtkSmartPointer<vtkPolyDataReader>::New();
    surfaceReader->SetFileName(inputSurfaceFile.c_str() );
    surfaceReader->Update();

    surface = surfaceReader->GetOutput();
    }

  // give the number of genus if demanded
  if( genusNumber )
    {
    int iNumberOfConnectedComponents = 1;
    int iNumberOfPoints = surface->GetPoints()->GetNumberOfPoints();
    int iNumberOfCells =  surface->GetNumberOfCells();

    vtkExtractEdges *extractEdges = vtkExtractEdges::New();
    extractEdges->SetInput( surface );
    extractEdges->Update();

    int iNumberOfEdges = extractEdges->GetOutput()->GetNumberOfLines();
    int iEulerCharacteristic = iNumberOfPoints - iNumberOfEdges + iNumberOfCells;

    std::cout << "number of connected components = " << iNumberOfConnectedComponents << std::endl;
    std::cout << "number of points = " << iNumberOfPoints << std::endl;
    std::cout << "number of triangles = " << iNumberOfCells << std::endl;
    std::cout << "number of edges = " << iNumberOfEdges << std::endl;
    std::cout << "Euler characteristic (sphere==2) = " << iEulerCharacteristic << std::endl;
    }

  // convert vtkPolyData to itkQuadEdgeMesh
  typedef float MeshPixelType;
  const unsigned int Dimension = 3;

  typedef itk::QuadEdgeMesh<MeshPixelType, Dimension> MeshType;

  typedef itk::ConvertVTKToQuadEdgeMeshFilter<MeshType> convertFilterType;
  convertFilterType::Pointer convertor = convertFilterType::New();
  convertor->SetPolyData(surface);
  convertor->Update();

  MeshType::Pointer mesh = convertor->GetOutput();
  mesh->DisconnectPipeline();

  // decimate the surface if requires
  if( decimateSurface )
    {
    std::cout << "Decimate the surface." << std::endl;

    typedef itk::NumberOfFacesCriterion<MeshType> CriterionType;
    CriterionType::Pointer criterion = CriterionType::New();
    criterion->SetTopologicalChange( false );
    criterion->SetNumberOfElements( numberOfElements );

    typedef itk::SquaredEdgeLengthDecimationQuadEdgeMeshFilter<MeshType, MeshType,
                                                               CriterionType> DecimationType;
    DecimationType::Pointer decimator = DecimationType::New();
    decimator->SetInput( convertor->GetOutput() );
    decimator->SetCriterion( criterion );
    decimator->Update();
    mesh = decimator->GetOutput();
    mesh->DisconnectPipeline();
    }

  // smooth the surface if requires
  if( smoothSurface )
    {
    std::cout << "Smooth the surface." << std::endl;

    itk::OnesMatrixCoefficients<MeshType> coeff0;

    typedef itk::SmoothingQuadEdgeMeshFilter<MeshType, MeshType> SmoothingType;
    SmoothingType::Pointer smoother = SmoothingType::New();
    smoother->SetInput( mesh );
    smoother->SetNumberOfIterations( numIterations );
    smoother->SetRelaxationFactor( relaxationFactor );
    smoother->SetCoefficientsMethod( &coeff0 );
    smoother->Update();

    mesh = smoother->GetOutput();
    mesh->DisconnectPipeline();
    }

  // write the output surface
  typedef itk::VTKPolyDataWriter<MeshType> WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetInput( mesh );
  writer->SetFileName( outputSurface.c_str() );
  writer->Update();

  return EXIT_SUCCESS;
}
