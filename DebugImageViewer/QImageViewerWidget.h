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
#ifndef QImageViewerWidget_h
#define QImageViewerWidget_h
#include <QWidget>
#include <vtkRenderer.h>

#include "QVTKWidget.h"
#include <list>

class vtkImageViewer2;

class QImageViewerWidget : public QVTKWidget
{
  Q_OBJECT
public slots:
  void SetSlice(int slice);

  void SetSliceOrientationToXY();

  void SetSliceOrientationToXZ();

  void SetSliceOrientationToYZ();

public:
  QImageViewerWidget(QWidget *parent = NULL);
  ~QImageViewerWidget();

  double GetColorWindow();

  double GetColorLevel();

  vtkRenderer * GetRenderer();

  void SetColorWindow(double s);

  void SetColorLevel(double s);

  void SetInput(vtkImageData *in);

  int * GetSliceRange();

  void Render();

  int GetSlice();

  int GetSliceMin();

  int GetSliceMax();

private:
  vtkImageViewer2* m_ImageViewer;
};

#endif // QImageViewerWidget_h
