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
#include <Qt>
#include <QMouseEvent>
#include <QInputEvent>
#include "QImageViewerWidget.h"
#include "vtkCamera.h"
#include "vtkImageViewer2.h"
#include "vtkRenderer.h"
#include "vtkRenderWindow.h"
#include "vtkImageViewer2.h"
#include "vtkInteractorStyleImage.h"
#include "vtkCommand.h"

QImageViewerWidget::QImageViewerWidget(QWidget *parent) : QVTKWidget(parent),
  m_ImageViewer(0)
{
  this->m_ImageViewer = vtkImageViewer2::New();
  this->SetRenderWindow(this->m_ImageViewer->GetRenderWindow() );
  this->m_ImageViewer->SetupInteractor
    (this->m_ImageViewer->GetRenderWindow()->GetInteractor() );
}

QImageViewerWidget::
~QImageViewerWidget()
{
}

void
QImageViewerWidget::SetSlice(int slice)
{
  this->m_ImageViewer->SetSlice(slice);
}

void
QImageViewerWidget::Render()
{
  this->m_ImageViewer->Render();
}

double
QImageViewerWidget::GetColorWindow()
{
  return this->m_ImageViewer->GetColorWindow();
}

double
QImageViewerWidget::GetColorLevel()
{
  return this->m_ImageViewer->GetColorLevel();
}

vtkRenderer *
QImageViewerWidget::GetRenderer()
{
  return this->m_ImageViewer->GetRenderer();
}

void
QImageViewerWidget::SetColorWindow(double s)
{
  return this->m_ImageViewer->SetColorWindow(s);
}

void
QImageViewerWidget::SetColorLevel(double s)
{
  return this->m_ImageViewer->SetColorLevel(s);
}

void
QImageViewerWidget::SetInput(vtkImageData *in)
{
  return this->m_ImageViewer->SetInput(in);
}

void
QImageViewerWidget::SetSliceOrientationToXY()
{
  this->m_ImageViewer->SetSliceOrientationToXY();
  vtkCamera *cam = this->m_ImageViewer->GetRenderer()->GetActiveCamera();
  cam->SetViewUp(0.0, -1.0, 0.0);
  this->m_ImageViewer->Render();
}

void
QImageViewerWidget::SetSliceOrientationToXZ()
{
  this->m_ImageViewer->SetSliceOrientationToXZ();
}

void
QImageViewerWidget::SetSliceOrientationToYZ()
{
  this->m_ImageViewer->SetSliceOrientationToYZ();
}

int *
QImageViewerWidget::GetSliceRange()
{
  return this->m_ImageViewer->GetSliceRange();
}

int
QImageViewerWidget::GetSlice()
{
  return this->m_ImageViewer->GetSlice();
}

int
QImageViewerWidget::GetSliceMin()
{
  return this->m_ImageViewer->GetSliceMin();
}

int
QImageViewerWidget::GetSliceMax()
{
  return this->m_ImageViewer->GetSliceMax();
}
