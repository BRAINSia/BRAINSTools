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
/*
 * Author: Wei Lu
 * at Psychiatry Imaging Lab,
 * University of Iowa Health Care 2010
 */

#ifndef _QSliceViewer_H
#define _QSliceViewer_H

#include "vtkSphereSource.h"
#include "vtkPolyDataMapper2D.h"
#include "vtkProperty2D.h"
#include "vtkCamera.h"
#include "vtkActor2D.h"
#include "vtkActor2DCollection.h"
#include "vtkRenderer.h"
#include "vtkRendererCollection.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"

#include "QVTKWidget.h"
#include <QObject>
#include <QWheelEvent>
#include <QListWidget>

#include <iostream>

#include "itkMacro.h" //Needed for nullptr

class QSliceViewer : public QVTKWidget
{
  Q_OBJECT
public:
  QSliceViewer( int type, QWidget * myParent = nullptr )
    : QVTKWidget( myParent )
    , m_bound( nullptr )
  {
    m_actors = vtkActor2DCollection::New();
    m_actor = nullptr;
    m_type = type;
    m_cPos = -1.0;
    m_color = 0;
    m_r = 0;
  }

  void
  SetBound( double * bound )
  {
    m_bound = bound;
  }

  double *
  GetBound()
  {
    return m_bound;
  }

  void
  createLabel( double * labelPos );

public slots:

  void
  createLabelSlot();

  void
  switchLabelSlot();

  void
  moveLabelSlot( double * labelPos ); // labelPos is a ratio

  void
  deleteLabelSlot();

  void
  deleteAllLabelSlot();

  void
  deleteLabelMouseSlot( QListWidgetItem * item ); // a mouse version wrap

  // for deleteLabelSlot

  void
  pickLabelSlot( QListWidgetItem * item ); // deal with click signal from

  // list

  void
  wheelSlot( double * labelPos ); // handle wheeling event

  void
  visibilityUpdate( int * table ); // update labels according to their

  // visibility
signals:
protected:
  vtkActor2DCollection * m_actors;

  vtkActor2D * m_actor;

  // slice viewer type: axial, sagittal, or coronal
  int m_type;

  // Pointer to the physical bound ( = ordered extent )
  double * m_bound;

  // last camera postion, helping to account wheeling effect
  double m_cPos;

  // color seed
  int m_color;

  // internal ratio for screen size change
  double m_r;

private:
  void
  GenSphere(); // generate a sphere

  void
  Highlight(); // highlight the current actor
};

#endif
