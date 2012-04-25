/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkTestingMacros.h,v $
  Language:  C++
  Date:      $Date: 2009-05-09 17:40:20 $
  Version:   $Revision: 1.1 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/

#ifndef __RegistrationMonitor_h
#define __RegistrationMonitor_h

#include "itkCommand.h"

#include "vtkPolyData.h"
#include "vtkPolyDataMapper.h"
#include "vtkActor.h"
#include "vtkProperty.h"
#include "vtkRenderer.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkMatrix4x4.h"
#include "vtkSmartPointer.h"
#include "vtkWindowToImageFilter.h"
#include "vtkJPEGWriter.h"
#include "vtkTextProperty.h"
#include "vtkTextActor.h"

/** \class RegistrationMonitor
 *  This class provides a VTK visualization pipeline configured for monitoring
 *  the progress of a registration process.
 */
class RegistrationMonitor
{
public:

  typedef RegistrationMonitor Self;

  RegistrationMonitor();
  virtual ~RegistrationMonitor();

  void SetFixedSurface( vtkPolyData* surface );

  void SetMovingSurface( vtkPolyData* surface );

  void SetNumberOfIterationsPerUpdate( unsigned int number );

  void Observe( itk::Object * processObject );

  void SetVerbose( bool );

  void SetOverlapMeshes( bool );

  void SetScreenShotsBaseFileName( const char * screenShotFileName );

  void SetBaseAnnotationText( const char * text );

  void SetCameraAzimuthAngle( double );

  void SetCameraElevationAngle( double );

  void SetCameraZoomFactor( double );

  void SetScalarRange( double minimun, double maximum );

protected:

  virtual void RefreshRendering();

  virtual void UpdateDataBeforeRendering() = 0;

  virtual void SetMovingActorMatrix( vtkMatrix4x4 * matrix );

  virtual void SetFixedActorMatrix( vtkMatrix4x4 * matrix );

  vtkPoints * GetFixedSurfacePoints();

  vtkPoints * GetMovingSurfacePoints();

  virtual void MarkFixedSurfaceAsModified();

  virtual void MarkMovingSurfaceAsModified();

  virtual void PrintOutUpdateMessage();

  itk::Object * GetObservedObject();

private:

  // These methods will only be called by the Observer
  virtual void Update();

  virtual void StartVisualization();

  virtual void RenderAndSaveScreenShot();

  void RemovePreviousObservers();

  void UpdateAnnotationText();

  vtkSmartPointer<vtkPolyData>       FixedSurface;
  vtkSmartPointer<vtkActor>          FixedActor;
  vtkSmartPointer<vtkProperty>       FixedProperty;
  vtkSmartPointer<vtkPolyDataMapper> FixedMapper;

  vtkSmartPointer<vtkPolyData>       MovingSurface;
  vtkSmartPointer<vtkActor>          MovingActor;
  vtkSmartPointer<vtkProperty>       MovingProperty;
  vtkSmartPointer<vtkPolyDataMapper> MovingMapper;

  // Visualization pipeline
  vtkSmartPointer<vtkRenderer>               FixedRenderer;
  vtkSmartPointer<vtkRenderer>               MovingRenderer;
  vtkSmartPointer<vtkRenderWindow>           RenderWindow;
  vtkSmartPointer<vtkRenderWindowInteractor> RenderWindowInteractor;
  vtkSmartPointer<vtkWindowToImageFilter>    WindowToImageFilter;

  vtkSmartPointer<vtkTextProperty> TextProperty;
  vtkSmartPointer<vtkTextActor>    TextActor;

  vtkSmartPointer<vtkJPEGWriter> Writer;

  typedef itk::SimpleMemberCommand<Self> ObserverType;

  ObserverType::Pointer IterationObserver;
  ObserverType::Pointer StartObserver;

  typedef std::vector<unsigned long> ObserverTagsArrayType;

  ObserverTagsArrayType ObserverTags;

  unsigned int CurrentIterationNumber;
  unsigned int NumberOfIterationsPerUpdate;

  bool Verbose;

  bool OverlapMeshes;

  itk::Object::Pointer ObservedObject;

  std::string ScreenShotsBaseFileName;

  std::string BaseAnnotationText;

  double CameraElevationAngle;
  double CameraAzimuthAngle;
  double CameraZoomFactor;

  double RangeMinimum;
  double RangeMaximum;
};

#endif
