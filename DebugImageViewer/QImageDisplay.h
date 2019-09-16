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
#ifndef QImageDisplay_h
#define QImageDisplay_h
#include "itkImage.h"
#include <QWidget>
#include <QString>
#include <string>

class QSlider;
class QImageViewerWidget;
class vtkKWImage;

class QImageDisplay : public QWidget
{
  Q_OBJECT
public slots:
  void
  viewChanged(const QString & newView);

public:
  using ReadImageType = itk::Image<float, 3>;
  using ImageType = itk::Image<unsigned char, 3>;
  using ImagePointer = ImageType::Pointer;
  QImageDisplay(QWidget * parent = 0);
  void
  SetImage(const std::string & fileName);

  void
  SetImage(ImagePointer & image);

  void
  SetBlankImage();

private:
  void
  SetSliceScaleRange();

  QImageViewerWidget * m_ImageViewer;
  vtkKWImage *         m_Image;
  QSlider *            m_Slider;
};

#endif // QImageDisplay_h
