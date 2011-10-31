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
  void viewChanged(const QString & newView);

public:
  typedef itk::Image<float, 3>         ReadImageType;
  typedef itk::Image<unsigned char, 3> ImageType;
  typedef ImageType::Pointer           ImagePointer;
  QImageDisplay(QWidget *parent = 0);
  void SetImage(const std::string & fileName);

  void SetImage(ImagePointer & image);

  void SetBlankImage();

private:
  void SetSliceScaleRange();

  QImageViewerWidget* m_ImageViewer;
  vtkKWImage*         m_Image;
  QSlider*            m_Slider;
};

#endif // QImageDisplay_h
