/*=auto=========================================================================

Portions (c) Copyright 2006 Brigham and Women's Hospital (BWH) All Rights Reserved.

See Doc/copyright/copyright.txt
or http://www.slicer.org/copyright/copyright.txt for details.

Program:   3D Slicer
Module:    $RCSfile: GenusZeroImageFilter.cxx,v $
Date:      $Date: 2009/08/17 16:12:02 $
Version:   $Revision: 1.4 $

=========================================================================auto=*/

#include "GenusZeroImageFilterCLP.h"
#include "vtkImageGenus0MarchingCubes.h"
#include "itkImageToVTKImageFilter.h"
#include "itkVTKImageToImageFilter.h"

#include "itkLabelExtracterImageFilter.h"

#include <itkImage.h>
#include <itkScalarConnectedComponentImageFilter.h>

#include "vtkImageGenus0MarchingCubes.h"

#include <vtkITKArchetypeImageSeriesReader.h>
#include <vtkITKArchetypeImageSeriesScalarReader.h>
#include <vtkImageData.h>
#include <vtkPointData.h>
#include <vtkTransform.h>
#include <vtkTransformPolyDataFilter.h>
#include <vtkPolyDataWriter.h>
#include <vtkImageChangeInformation.h>

#include <ModuleDescriptionParser.h>
#include <ModuleDescription.h>
#include <ModuleParameterGroup.h>
#include <ModuleParameter.h>

#include <vtkExtractEdges.h>

#include <vtkPolyDataWriter.h>
#include <vtkSTLWriter.h>
#include <vtkImageCast.h>

#include <vtkITKImageWriter.h>

#include <map>
#include <string>

int
main( int argc, char * argv[] )
{
  PARSE_ARGS;

  bool debug = true;

  // vtk and helper variables
  vtkITKArchetypeImageSeriesReader * reader = nullptr;
  vtkImageData *                     image;
  vtkTransform *                     transformIJKtoRAS = nullptr;
  vtkTransformPolyDataFilter *       transformer = nullptr;
  vtkPolyDataWriter *                writer = nullptr;

  using ImageType = itk::Image< unsigned short, 3 >;
  using VTKToITKImageFilterType = itk::VTKImageToImageFilter< ImageType >;
  using ITKToVTKImageFilterType = itk::ImageToVTKImageFilter< ImageType >;
  using ConnectedComponentImageFilterType = itk::ScalarConnectedComponentImageFilter< ImageType, ImageType >;
  using LabelExtracterImageFilterType = itk::LabelExtracterImageFilter< ImageType, ImageType >;
  // check for the input file
  FILE * infile;
  infile = fopen( inputVolume.c_str(), "r" );
  if ( infile == nullptr )
  {
    std::cerr << "ERROR: cannot open input volume file " << inputVolume << endl;
    return EXIT_FAILURE;
  }
  fclose( infile );

  // Read the file
  reader = vtkITKArchetypeImageSeriesScalarReader::New();
  reader->SetArchetype( inputVolume.c_str() );
  reader->SetOutputScalarTypeToDouble();
  reader->SetDesiredCoordinateOrientationToNative();
  reader->SetUseNativeOriginOn();
  reader->Update();

  std::cout << "Done reading the file " << inputVolume << endl;

  // change the image information to a canonical form that we can easily process
  // in the vtk filter

  vtkImageChangeInformation * ici = vtkImageChangeInformation::New();
  ici->SetInputData( reader->GetOutput() );
  ici->SetOutputSpacing( 1, 1, 1 );
  ici->SetOutputOrigin( 0, 0, 0 );
  ici->Update();

  image = ici->GetOutput();
  // image->UpdateInformation();

  // Get the dimensions, marching cubes only works on 3d
  int extents[6];
  image->GetExtent( extents );
  if ( debug )
  {
    std::cout << "Image data extents: " << extents[0] << " " << extents[1] << " " << extents[2] << " " << extents[3]
              << " " << extents[4] << " " << extents[5] << endl;
  }
  if ( extents[0] == extents[1] || extents[2] == extents[3] || extents[4] == extents[5] )
  {
    std::cerr << "The volume is not 3D.\n";
    std::cerr << "\tImage data extents: " << extents[0] << " " << extents[1] << " " << extents[2] << " " << extents[3]
              << " " << extents[4] << " " << extents[5] << endl;
    return EXIT_FAILURE;
  }

  // Get the RAS to IJK matrix and invert it to get the IJK to RAS which will
  // need
  // to be applied to the model as it will be built in pixel space

  transformIJKtoRAS = vtkTransform::New();
  transformIJKtoRAS->SetMatrix( reader->GetRasToIjkMatrix() );
  if ( debug )
  {
    std::cout << "RasToIjk matrix from file = ";
    transformIJKtoRAS->GetMatrix()->Print( std::cout );
  }
  transformIJKtoRAS->Inverse();

  // set up the genus zero marching cubes filter

  vtkImageGenus0MarchingCubes * marchingcubes = vtkImageGenus0MarchingCubes::New();
  /*    vtkPluginFilterWatcher watchMCubes(marchingcubes,
                                     "Marching Cubes",
                                     CLPProcessInformation,
                                     1.0/7.0, 0.0);*/

  // pass the command line arguments to the filter

  marchingcubes->SetInputData( ici->GetOutput() );

  // check the supported cases and abort in case of an unsupported
  // case

  if ( !( connectivityModel == 18 && computeSurface ) && !( connectivityModel == 6 ) )
  {
    std::cerr << "Mode: "
              << " connectivityModel =  " << connectivityModel << "; computeSurface = " << computeSurface
              << "  not supported in current implementation. ABORT." << std::endl;
    return EXIT_FAILURE;
  }

  std::cout << "biggestComponent = " << biggestComponent << std::endl;
  std::cout << "connectivityModel = " << connectivityModel << std::endl;
  std::cout << "cutLoops = " << cutLoops << std::endl;
  std::cout << "computeSurface = " << computeSurface << std::endl;

  if ( computeSurface )
  {
    marchingcubes->ComputeSurfaceOn();
  }
  else
  {
    marchingcubes->ComputeSurfaceOff();
  }

  if ( !biggestComponent )
  {
    marchingcubes->BiggestComponentOff();
  }
  else
  {
    marchingcubes->BiggestComponentOn();
  }

  if ( connectivityModel == 18 )
  {
    marchingcubes->Use18Connectivity();
  }
  else
  {
    marchingcubes->Use6Connectivity();
  }

  if ( connectedComponent )
  {
    marchingcubes->ConnectedComponentOn();
  }
  else
  {
    marchingcubes->ConnectedComponentOff();
  }

  if ( cutLoops )
  {
    marchingcubes->CutLoopsOn();
  }
  else
  {
    marchingcubes->CutLoopsOff();
  }

  marchingcubes->Update();

  if ( 1 ) // computeSurface )
  {

    // if the surface was computed, compute it's topology

    int iNumberOfConnectedComponents = marchingcubes->GetNumberOfConnectedComponents();

    // determine the Euler characteristic of this

    int iNumberOfPoints = marchingcubes->GetOutput()->GetPoints()->GetNumberOfPoints();
    int iNumberOfCells = marchingcubes->GetOutput()->GetNumberOfCells();

    // determine number of edges

    vtkExtractEdges * extractEdges = vtkExtractEdges::New();
    extractEdges->SetInputData( marchingcubes->GetOutput() );
    extractEdges->Update();

    int iNumberOfEdges = extractEdges->GetOutput()->GetNumberOfLines();
    int iEulerCharacteristic = iNumberOfPoints - iNumberOfEdges + iNumberOfCells;

    std::cout << "number of connected components = " << iNumberOfConnectedComponents << std::endl;
    std::cout << "number of points = " << iNumberOfPoints << std::endl;
    std::cout << "number of triangles = " << iNumberOfCells << std::endl;
    std::cout << "number of edges = " << iNumberOfEdges << std::endl;
    std::cout << "Euler characteristic (sphere==2) = " << iEulerCharacteristic << std::endl;

    if ( iEulerCharacteristic != 2 * iNumberOfConnectedComponents )
    {
      std::cerr << "Object deviates topologically from a sphere. Results may not be the desired one ..." << std::endl;
    }
  }

  // transform the vtk polydata back to RAS

  transformer = vtkTransformPolyDataFilter::New();
  /*    vtkPluginFilterWatcher watchTranformer(transformer,
                                         "Transformer",
                                         CLPProcessInformation,
                                         1.0/7.0, 4.0/7.0);*/
  transformer->SetInputData( marchingcubes->GetOutput() );

  transformer->SetTransform( transformIJKtoRAS );
  if ( debug )
  {
    std::cout << "Transforming using inversed matrix:\n";
    transformIJKtoRAS->GetMatrix()->Print( std::cout );
  }
  transformer->Update();

  // INFO: add progress
  // (transformer->GetOutput())->ReleaseDataFlagOn();

  // but for now we're just going to write it out

  writer = vtkPolyDataWriter::New();
  writer->SetInputData( transformer->GetOutput() );

  if ( computeSurface )
  {
    writer->SetFileName( vtkOutput.c_str() );
    // INFO: add progress
    writer->Write();
  }

  // now write out the topology corrected volume

  vtkITKImageWriter * imageWriter = vtkITKImageWriter::New();
  imageWriter->SetRasToIJKMatrix( reader->GetRasToIjkMatrix() );

  if ( !computeSurface )
  {

    vtkImageCast * imageCast = vtkImageCast::New();
    imageCast->SetInputData( marchingcubes->GetCorrectedImageData() );
    imageCast->SetOutputScalarTypeToUnsignedShort();
    imageCast->Update();

    if ( extractFinalConnectedComponent )
    {
      std::cout << "Extracting the largest connected component of the resulting image volume." << std::endl;

      // need to go back and forth from vtk to itk (and back)

      VTKToITKImageFilterType::Pointer           vtk2itkImageFilter = VTKToITKImageFilterType::New();
      ITKToVTKImageFilterType::Pointer           itk2vtkImageFilter = ITKToVTKImageFilterType::New();
      ConnectedComponentImageFilterType::Pointer connectedComponentFilter = ConnectedComponentImageFilterType::New();
      LabelExtracterImageFilterType::Pointer     labelExtracter = LabelExtracterImageFilterType::New();

      vtk2itkImageFilter->SetInputData( imageCast->GetOutput() );
      vtk2itkImageFilter->Update();

      // get the connected components

      connectedComponentFilter->SetDistanceThreshold( 0 );
      connectedComponentFilter->FullyConnectedOff();

      if ( connectivityModel != 6 )
      {
        std::cerr << "WARNING: Only 6 connectivity is supported for the extraction of the largest connected component "
                     "of the result. DEFAULTING to 6 connectivity."
                  << std::endl;
      }

      connectedComponentFilter->SetInput( vtk2itkImageFilter->GetOutput() );
      connectedComponentFilter->Update();

      // now find the largest connected component and extract it

      // iterate throught the output and figure out how many elements we have,
      // do this with a map

      std::map< int, int > ccsMap;

      using ConstIteratorType = itk::ImageRegionConstIterator< ImageType >;
      ImageType::ConstPointer ccs = connectedComponentFilter->GetOutput();

      ConstIteratorType it( ccs, ccs->GetLargestPossibleRegion() );
      for ( it.GoToBegin(); !it.IsAtEnd(); ++it )
      {

        int ip = (int)round( it.Get() );

        if ( ccsMap.count( ip ) == 0 )
        {
          ccsMap[ip] = 1;
        }
        else
        {
          ccsMap[ip] += 1;
        }
      }

      // determine the largest connected component

      // the first one (1) will be the background

      int iLargestComponentId = -1;
      int iNumber = 0;

      std::cout << "Connected components found = " << std::endl;

      std::map< int, int >::iterator iterMap;
      for ( iterMap = ccsMap.begin(); iterMap != ccsMap.end(); iterMap++ )
      {
        int iId = ( *iterMap ).first;
        int iCount = ( *iterMap ).second;

        std::cout << iId << " -> " << iCount << std::endl;

        if ( iId > 1 && iCount > iNumber )
        {
          iNumber = iCount;
          iLargestComponentId = iId;
        }
      }

      std::cout << "Extracting component " << iLargestComponentId << std::endl;

      // creating the change map

      LabelExtracterImageFilterType::ChangeMapType changeMap;

      changeMap[iLargestComponentId] = 1;

      labelExtracter->SetChangeMap( changeMap );
      labelExtracter->SetInput( connectedComponentFilter->GetOutput() );
      labelExtracter->Update();

      itk2vtkImageFilter->SetInputData( labelExtracter->GetOutput() );
      // itk2vtkImageFilter->Update();

      // and write out the result

      imageWriter->SetInputData( itk2vtkImageFilter->GetOutput() );

      imageWriter->SetFileName( outputVolume.c_str() );
#if ITK_VERSION_MAJOR >= 5
      imageWriter->SetUseCompression( 1 );
#endif
      imageWriter->Write();
    }
    else // just write out everything there is
    {

      imageWriter->SetInputData( imageCast->GetOutput() );

      imageWriter->SetFileName( outputVolume.c_str() );
#if ITK_VERSION_MAJOR >= 5
      imageWriter->SetUseCompression( 1 );
#endif
      imageWriter->Write();
    }
  }

  // Cleanup
  if ( reader )
  {
    reader->Delete();
  }
  if ( ici )
  {
    ici->Delete();
  }
  if ( transformIJKtoRAS )
  {
    transformIJKtoRAS->Delete();
  }
  if ( marchingcubes )
  {
    marchingcubes->Delete();
  }
  if ( transformer )
  {
    transformer->Delete();
  }
  if ( writer )
  {
    writer->Delete();
  }
  if ( imageWriter )
  {
    imageWriter->Delete();
  }
  return EXIT_SUCCESS;
}
