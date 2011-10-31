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
