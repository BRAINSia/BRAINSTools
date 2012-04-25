/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkFilterWatcher.h,v $
  Language:  C++
  Date:      $Date: 2007-01-29 14:42:11 $
  Version:   $Revision: 1.15 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/

#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif

#include "vtkRenderer.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkActor.h"
#include "vtkCamera.h"
#include "vtkPolyData.h"
#include "vtkPolyDataMapper.h"
#include "vtkLookupTable.h"
#include "vtkColorTransferFunction.h"
#include "vtkContourFilter.h"
#include "vtkPolyDataReader.h"
#include "vtkProperty.h"
#include "vtkDataSet.h"
#include "vtkWindowToImageFilter.h"
#include "vtkPNGWriter.h"
#include "vtkScalarBarActor.h"
#include "vtkInteractorStyleTrackballCamera.h"

#include "vtksys/SystemTools.hxx"

int main( int argc, char* argv[] )
{
  if( argc < 2 )
    {
    std::cerr << "Usage: " << std::endl;
    std::cerr << argv[0] << " InputMesh ";
    std::cerr << " [OutputPNGScreenshot]" << std::endl;
    return EXIT_FAILURE;
    }

  vtkPolyDataReader * surfaceReader = vtkPolyDataReader::New();
  surfaceReader->SetFileName( argv[1] );
  surfaceReader->Update();

  vtkPolyData * surface = surfaceReader->GetOutput();

  // Create a renderer, render window, and render window interactor to
  // display the results.
  vtkRenderer *               renderer = vtkRenderer::New();
  vtkRenderWindow *           renWin = vtkRenderWindow::New();
  vtkRenderWindowInteractor * iren = vtkRenderWindowInteractor::New();

  renWin->SetSize( 512, 512 );
  renWin->AddRenderer(renderer);
  iren->SetRenderWindow(renWin);

  renderer->SetBackground( 1.0, 1.0, 1.0 );

  vtkPolyDataMapper * polyMapper = vtkPolyDataMapper::New();
  vtkActor *          polyActor  = vtkActor::New();

  polyActor->SetMapper( polyMapper );
  polyMapper->SetInput( surface );
  polyMapper->ScalarVisibilityOn();
  polyMapper->SetScalarModeToUsePointData();
  polyMapper->SetColorModeToMapScalars();

  vtkLookupTable * lookUpTable = vtkLookupTable::New();
  lookUpTable->SetRampToLinear();
  lookUpTable->SetNumberOfColors( 256 );
  lookUpTable->SetHueRange(0.667, 0.0);
  lookUpTable->SetRange(0.0, 1.0);
  lookUpTable->Build();

  vtkColorTransferFunction * colorFunction = vtkColorTransferFunction::New();
  colorFunction->SetColorSpaceToHSV();
  colorFunction->AddRGBPoint( 0.0, 0.0, 0.0, 1.0 );
  colorFunction->AddRGBPoint( 1.0, 1.0, 0.0, 0.0 );

/*  colorFunction->AddHSVPoint( 0.0, 0.667, 1.0, 1.0);
  colorFunction->AddHSVPoint( 0.5, 0.333, 1.0, 0.0);
  colorFunction->AddHSVPoint( 1.0, 0.000, 1.0, 1.0);
*/

  polyMapper->SetLookupTable( lookUpTable );
  // polyMapper->SetLookupTable( colorFunction );
  polyMapper->SetScalarRange( 0.0, 1.0 );

  // vtkScalarBarActor * scalarBarActor = vtkScalarBarActor::New();
  // scalarBarActor->SetTitle("Density");
  // scalarBarActor->SetLookupTable( lookUpTable );
  // renderer->AddActor2D(scalarBarActor);

  vtkProperty * property = vtkProperty::New();
  property->SetAmbient(0.1);
  property->SetDiffuse(0.8);
  property->SetSpecular(0.1);
  property->SetLineWidth(2.0);
  property->SetRepresentationToSurface();

  polyActor->SetProperty( property );

  renderer->AddActor( polyActor );

  vtkInteractorStyleTrackballCamera * interactorStyle = vtkInteractorStyleTrackballCamera::New();
  iren->SetInteractorStyle( interactorStyle );

  double bounds[6];
  surface->GetBounds( bounds );

  double centerOfPolydata[3];

  centerOfPolydata[0] = ( bounds[0] + bounds[1] ) / 2.0;
  centerOfPolydata[1] = ( bounds[2] + bounds[3] ) / 2.0;
  centerOfPolydata[2] = ( bounds[4] + bounds[5] ) / 2.0;

  double cameraPosition[3];

  double zoomFactor = 1.5;

  if( argc <= 2 )
    {
    // Bring up the render window and begin interaction.
    vtkCamera * camera = renderer->GetActiveCamera();

    cameraPosition[0] = centerOfPolydata[0] + 200.0;
    cameraPosition[1] = centerOfPolydata[1] - 100.0;
    cameraPosition[2] = centerOfPolydata[2] + 500.0;

    camera->SetPosition( cameraPosition );
    camera->SetFocalPoint( centerOfPolydata );
    camera->SetViewUp( 0, 1, 0 );
    camera->SetParallelProjection( true );
    renderer->ResetCamera();
    camera->Zoom( zoomFactor );
    renWin->Render();
    iren->Start();
    }

  if( argc > 2 )
    {
    vtkWindowToImageFilter * windowToImageFilter = vtkWindowToImageFilter::New();

    vtkPNGWriter * screenShotWriter = vtkPNGWriter::New();

    vtkCamera * camera = renderer->GetActiveCamera();

    cameraPosition[0] = centerOfPolydata[0] + 200.0;
    cameraPosition[1] = centerOfPolydata[1] - 100.0;
    cameraPosition[2] = centerOfPolydata[2] + 500.0;

    camera->SetPosition( cameraPosition );
    camera->SetFocalPoint( centerOfPolydata );
    camera->SetViewUp( 0, 1, 0 );
    camera->SetParallelProjection( true );
    renderer->ResetCamera();
    camera->Zoom( zoomFactor );

    windowToImageFilter->SetInput( renWin );
    windowToImageFilter->Update();

    screenShotWriter->SetInput( windowToImageFilter->GetOutput() );
    screenShotWriter->SetFileName( argv[2] );

    renWin->Render();
    screenShotWriter->Write();

    screenShotWriter->SetInput( NULL );
    windowToImageFilter->SetInput( NULL );

    windowToImageFilter->Delete();
    screenShotWriter->Delete();
    }

  // Release all VTK components
  polyActor->Delete();
  property->Delete();
  // lookUpTable->Delete();
  // colorFunction->Delete();
  // scalarBarActor->Delete();
  polyMapper->Delete();
  interactorStyle->Delete();
  renWin->Delete();
  renderer->Delete();
  iren->Delete();

  // Destroy the objects
  surfaceReader->Delete();

  return EXIT_SUCCESS;
}
