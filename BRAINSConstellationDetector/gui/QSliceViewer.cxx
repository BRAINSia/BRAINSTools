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
#include "vtkVersion.h"
#include "QSliceViewer.h"

#include <iostream>

void
QSliceViewer::createLabel(double * labelPos)
{
  createLabelSlot();

  double * bound = m_bound;
  double * labelPosRatio = new double[3];
  labelPosRatio[0] = (labelPos[0] - bound[0]) / (bound[1] - bound[0]);
  labelPosRatio[1] = (labelPos[1] - bound[2]) / (bound[3] - bound[2]);
  labelPosRatio[2] = (labelPos[2] - bound[4]) / (bound[5] - bound[4]);

  double R; // physical width/height ratio
  if (m_type == 2)
  { // sagittal view mode ( ASL )
    R = (bound[3] - bound[2]) / (bound[5] - bound[4]);
  }
  else if (m_type == 3)
  { // coronal view mode ( LSA )
    R = (bound[1] - bound[0]) / (bound[5] - bound[4]);
  }
  else
  { // axial view mode ( LAS )
    R = (bound[1] - bound[0]) / (bound[3] - bound[2]);
  }

  // by observation, linear approximation
  double estimateCoef = 1. / sqrt(1.0 + R * R);
  labelPosRatio[0] = estimateCoef * (labelPosRatio[0] - .5) + .5;
  labelPosRatio[1] = estimateCoef * (labelPosRatio[1] - .5) + .5;
  labelPosRatio[2] = estimateCoef * (labelPosRatio[2] - .5) + .5;

  moveLabelSlot(labelPosRatio);

  // hide all actors
  m_actors->InitTraversal();
  vtkActor2D * actor = m_actors->GetNextActor2D();
  while (actor != nullptr)
  {
    actor->GetProperty()->SetOpacity(0.0);
    actor = m_actors->GetNextActor2D();
  }
}

void
QSliceViewer::createLabelSlot()
{
  GenSphere();
  m_actors->AddItem(m_actor);
  Highlight();
}

void
QSliceViewer::deleteLabelSlot()
{
  if (m_actor != nullptr)
  {
    this->GetRenderWindow()->GetRenderers()->GetFirstRenderer()->RemoveActor(m_actor);
    m_actors->RemoveItem(m_actor);
    m_actor->Delete();
    this->GetRenderWindow()->Render();
    m_actors->InitTraversal();
    m_actor = m_actors->GetNextActor2D();
    Highlight();
  }
}

void
QSliceViewer::deleteAllLabelSlot()
{
  m_actors->InitTraversal();
  m_actor = m_actors->GetNextActor2D();
  while (m_actor != nullptr)
  {
    this->GetRenderWindow()->GetRenderers()->GetFirstRenderer()->RemoveActor(m_actor);
    m_actors->RemoveItem(m_actor);
    m_actor->Delete();
    m_actor = m_actors->GetNextActor2D();
  }

  // reset the color seed
  m_color = 0;
  this->GetRenderWindow()->Render();
}

void
QSliceViewer::deleteLabelMouseSlot(QListWidgetItem *)
{
  deleteLabelSlot();
}

void
QSliceViewer::switchLabelSlot()
{
  if (m_actor != nullptr)
  {
    m_actor = m_actors->GetNextActor2D();
    if (m_actor == nullptr)
    {
      m_actors->InitTraversal();            // move current list pointer to the
                                            // beginning
      m_actor = m_actors->GetNextActor2D(); // set current label to the first
                                            // actor
    }
    Highlight();
  }
}

void
QSliceViewer::moveLabelSlot(double * labelPos)
{
  if (m_actor != nullptr) // we only need to update current label
  {
    // get current position for each viewer
    int              currPos[2];
    int              eventSize[2];
    QVTKInteractor * interactor = this->GetInteractor();
    interactor->GetSize(eventSize);
    double * bound = m_bound;
    double   cPos[3]; // camera position
    interactor->GetRenderWindow()->GetRenderers()->GetFirstRenderer()->GetActiveCamera()->GetPosition(cPos);

    // compute the screen size of image slice
    double R; // physical width/height ratio
    if (m_type == 2)
    { // sagittal view mode ( ASL )
      R = (bound[3] - bound[2]) / (bound[5] - bound[4]);
    }
    else if (m_type == 3)
    { // coronal view mode ( LSA )
      R = (bound[1] - bound[0]) / (bound[5] - bound[4]);
    }
    else
    { // axial view mode ( LAS )
      R = (bound[1] - bound[0]) / (bound[3] - bound[2]);
    }

    // by observation
    int h = 0.96 * eventSize[1] / sqrt(1.0 + R * R); // sreen height of image
                                                     // slice
    int w = R * h;                                   // screen width of image
                                                     // slice

    // wheeling effect
    if (m_cPos < 0.0)
    {
      m_cPos = cPos[2];
    }
    else if (m_cPos != cPos[2])
    {
      h *= m_cPos / cPos[2];
      w *= m_cPos / cPos[2];
    }

    // locate the label position
    if (m_type == 2)
    { // sagittal view mode ( ASL )
      currPos[0] = w * labelPos[1] - w / 2.0;
      currPos[1] = h * labelPos[2] - h / 2.0;
    }
    else if (m_type == 3)
    { // coronal view mode ( LSA )
      currPos[0] = w * labelPos[0] - w / 2.0;
      currPos[1] = h * labelPos[2] - h / 2.0;
    }
    else
    { // axial view mode ( LAS )
      currPos[0] = w * labelPos[0] - w / 2.0;
      currPos[1] = h * (1 - labelPos[1]) - h / 2.0;
    }

    // move sphere
    m_actor->SetPosition(currPos[0], currPos[1]);
    this->GetRenderWindow()->Render();
  }
}

void
QSliceViewer::wheelSlot(double * labelPos)
{
  if (m_actor != nullptr)
  {
    vtkActor2D * currActor = m_actor; // a backup
    m_actors->InitTraversal();
    vtkActor2D * actor = m_actors->GetNextActor2D();
    int          index = 0;        // index of actor
    double       labelPosRatio[3]; // label postion ratio for inter-viewer
                                   // communication
    double * bound = m_bound;

    while (actor != nullptr)
    {
      m_actor = actor;
      labelPosRatio[0] = (labelPos[3 * index] - bound[0]) / (bound[1] - bound[0]);
      labelPosRatio[1] = (labelPos[3 * index + 1] - bound[2]) / (bound[3] - bound[2]);
      labelPosRatio[2] = (labelPos[3 * index + 2] - bound[4]) / (bound[5] - bound[4]);
      moveLabelSlot(labelPosRatio);
      actor = m_actors->GetNextActor2D();
      ++index;
    }

    // restore list pointer position
    m_actor = currActor;
    m_actors->InitTraversal();
    actor = m_actors->GetNextActor2D();
    while (actor != m_actor)
    {
      actor = m_actors->GetNextActor2D();
    }

    this->GetRenderWindow()->Render();
  }
}

void
QSliceViewer::GenSphere()
{
  // create sphere geometry
  vtkSphereSource * sphere = vtkSphereSource::New();

  sphere->SetRadius(3.0);
  sphere->SetThetaResolution(20);
  sphere->SetPhiResolution(10);
  QVTKInteractor * interactor = this->GetInteractor();
  int              eventSize[2];
  interactor->GetSize(eventSize);
  double sphereCenter[3] = { eventSize[0] / 2.0, eventSize[1] / 2.0, 0 };
  sphere->SetCenter(sphereCenter);

  // map to graphics library
  vtkPolyDataMapper2D * map = vtkPolyDataMapper2D::New();
#if (VTK_MAJOR_VERSION < 6)
  map->SetInput(sphere->GetOutput());
#else
  map->SetInputData(sphere->GetOutput());
#endif

  // actor coordinates geometry, properties, transformation
  vtkActor2D * actor = vtkActor2D::New();
  actor->SetMapper(map);
  this->GetRenderWindow()->GetRenderers()->GetFirstRenderer()->AddActor(actor);

  m_actor = actor; // set current label to new created actor

  m_actor->GetProperty()->SetColor(m_color % 3 / 2 / 2.0 + (m_color * m_color) % 3 / 4.0,
                                   (m_color + 1) % 3 / 2 / 2.0 + m_color % 7 / 12.0,
                                   (m_color + 2) % 3 / 2 / 2.0 + m_color % 4 / 6.0);
  ++m_color;
}

void
QSliceViewer::Highlight()
{
  if (m_actor != nullptr) // we only need to update current label
  {
    // create sphere geometry
    vtkSphereSource * sphere1 = vtkSphereSource::New();
    sphere1->SetRadius(5.0);
    sphere1->SetThetaResolution(4);
    sphere1->SetPhiResolution(4);

    vtkSphereSource * sphere2 = vtkSphereSource::New();
    sphere2->SetRadius(3.0);
    sphere2->SetThetaResolution(20);
    sphere2->SetPhiResolution(10);

    QVTKInteractor * interactor = this->GetInteractor();
    int              eventSize[2];
    interactor->GetSize(eventSize);
    double sphereCenter[3] = { eventSize[0] / 2.0, eventSize[1] / 2.0, 0 };
    sphere1->SetCenter(sphereCenter);
    sphere2->SetCenter(sphereCenter);

    // map to graphics library
    vtkPolyDataMapper2D * map1 = vtkPolyDataMapper2D::New();
#if (VTK_MAJOR_VERSION < 6)
    map1->SetInput(sphere1->GetOutput());
#else
    map1->SetInputData(sphere1->GetOutput());
#endif
    vtkPolyDataMapper2D * map2 = vtkPolyDataMapper2D::New();
#if (VTK_MAJOR_VERSION < 6)
    map2->SetInput(sphere2->GetOutput());
#else
    map2->SetInputData(sphere2->GetOutput());
#endif

    m_actors->InitTraversal();
    vtkActor2D * actor = m_actors->GetNextActor2D();
    while (actor != nullptr)
    {
      double color[3];
      actor->GetProperty()->GetColor(color);
      double * myPos = actor->GetPosition();
      if (actor == m_actor)
      { // set size and shape for the current actor
        actor->SetMapper(map1);
      }
      else
      { // set any other actors to default size and shape
        actor->SetMapper(map2);
      }
      actor->GetProperty()->SetColor(color);
      actor->SetPosition(myPos[0], myPos[1]);
      actor = m_actors->GetNextActor2D();
    }

    // restore list pointer position
    m_actors->InitTraversal();
    actor = m_actors->GetNextActor2D();
    while (actor != m_actor)
    {
      actor = m_actors->GetNextActor2D();
    }

    this->GetRenderWindow()->Render();
  }
}

void
QSliceViewer::pickLabelSlot(QListWidgetItem * item)
{
  if (m_actor != nullptr) // we only need to update current label
  {
    m_actors->InitTraversal();

    int i;
    for (i = 0; i <= item->listWidget()->currentRow(); ++i)
    {
      m_actor = m_actors->GetNextActor2D();
    }
    m_actor->GetProperty()->SetOpacity(1.0);
    Highlight();
  }
}

void
QSliceViewer::visibilityUpdate(int * table)
{
  if (m_actor != nullptr)
  {
    m_actors->InitTraversal();
    vtkActor2D * actor = m_actors->GetNextActor2D();
    int          index = 0; // index of actor
    while (actor != nullptr)
    {
      actor = m_actors->GetNextActor2D();
      ++index;
    }

    int numActors = index;

    index = 0;
    m_actors->InitTraversal();
    actor = m_actors->GetNextActor2D();
    while (actor != nullptr)
    {
      if (table[numActors * ((m_type + 1) % 3) + index] == 0)
      {
        actor->GetProperty()->SetOpacity(0.0);
      }
      else
      {
        actor->GetProperty()->SetOpacity(1.0);
      }

      actor = m_actors->GetNextActor2D();
      ++index;
    }

    // restore list pointer position
    m_actors->InitTraversal();
    actor = m_actors->GetNextActor2D();
    while (actor != m_actor)
    {
      actor = m_actors->GetNextActor2D();
    }

    this->GetRenderWindow()->Render();
  }
}
