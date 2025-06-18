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
#include "QImageDisplay.h"
#include "QImageViewerWidget.h"
#include "itkIO.h"
#include "vtkKWImage.h"
#include "vtkImageData.h"
#include "DebugImageViewerClient.h"

#include <QGridLayout>
#include <QSlider>
#include <QComboBox>

void
QImageDisplay::SetSliceScaleRange()
{
  const int * range = this->m_ImageViewer->GetSliceRange();

  this->m_Slider->setRange(range[0], range[1]);
  int centerSlice = range[0] + ((range[1] - range[0]) / 2);
  this->m_Slider->setValue(centerSlice);
  this->m_ImageViewer->SetSlice(centerSlice);
}

void
QImageDisplay::viewChanged(const QString & newView)
{
  if (newView == "XY")
  {
    this->m_ImageViewer->SetSliceOrientationToXY();
  }
  else if (newView == "XZ")
  {
    this->m_ImageViewer->SetSliceOrientationToXZ();
  }
  else if (newView == "YZ")
  {
    this->m_ImageViewer->SetSliceOrientationToYZ();
  }
  this->SetSliceScaleRange();
}

void
QImageDisplay::SetImage(ImagePointer & image)
{
  if (this->m_Image != 0)
  {
    this->m_Image->Delete();
  }
  this->m_Image = vtkKWImage::New();
  this->m_Image->SetITKImageBase(image);
  vtkImageData * vtkImage = this->m_Image->GetVTKImage();
  double *       range = vtkImage->GetScalarRange();
  this->m_ImageViewer->SetInput(vtkImage);
  this->m_ImageViewer->SetColorWindow(range[1] - range[0]);
  this->m_ImageViewer->SetColorLevel(0.5 * (range[1] + range[0]));
  this->SetSliceScaleRange();
  vtkImage->Modified();
}

void
QImageDisplay::SetBlankImage()
{
  ImageType::SizeType    imageSize;
  ImageType::IndexType   imageIndex;
  ImageType::RegionType  imageRegion;
  ImageType::SpacingType imageSpacing;

  for (unsigned int i = 0; i < 3; i++)
  {
    imageSize[i] = 10;
    imageIndex[i] = 0;
    imageSpacing[i] = 1.0;
  }
  imageRegion.SetSize(imageSize);
  imageRegion.SetIndex(imageIndex);
  ImagePointer image = ImageType::New();
  image->SetSpacing(imageSpacing);
  image->SetRegions(imageRegion);
  image->Allocate();
  image->FillBuffer(0);
  this->SetImage(image);
}

void
QImageDisplay::SetImage(const std::string & fileName)
{
  ReadImageType::Pointer readImage =
    itkUtil::ReadImageAndOrient<ReadImageType>(fileName, itk::SpatialOrientationEnums::ITK_COORDINATE_ORIENTATION_RAI);
  ImageType::Pointer image = DebugImageViewerUtil::ScaleAndCast<ReadImageType, ImageType>(readImage, 0, 255);

  this->SetImage(image);
}

QImageDisplay::QImageDisplay(QWidget * parent)
  : QWidget(parent)
  , m_ImageViewer(0)
  , m_Image(0)
  , m_Slider(0)
{
  QGridLayout * layout = new QGridLayout(this);

  this->m_ImageViewer = new QImageViewerWidget(this);
  this->m_ImageViewer->setMinimumSize(275, 225);
  this->m_ImageViewer->SetSliceOrientationToXY();

  QComboBox * combo = new QComboBox(this);

  const char * viewNames[3] = { "XY", "XZ", "YZ" };
  for (unsigned i = 0; i < 3; i++)
  {
    combo->insertItem(i, tr(viewNames[i]));
  }

  connect(combo, SIGNAL(currentIndexChanged(const QString &)), this, SLOT(viewChanged(const QString &)));

  this->m_Slider = new QSlider(Qt::Vertical, this);
  this->m_Slider->setSingleStep(1);
  connect(this->m_Slider, SIGNAL(valueChanged(int)), this->m_ImageViewer, SLOT(SetSlice(int)));

  layout->addWidget(this->m_ImageViewer, 0, 0, 5, 5);
  layout->addWidget(combo, 5, 0, 6, 1);
  layout->addWidget(this->m_Slider, 0, 5, 5, 1);
  this->setLayout(layout);
}
