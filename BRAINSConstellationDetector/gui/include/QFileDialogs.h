/*
 * Author: Wei Lu
 * at Psychiatry Imaging Lab,
 * University of Iowa Health Care 2010
 */

#ifndef _QFileDialogs_H
#define _QFileDialogs_H

#include <QFileDialog>
#include <QObject>

class QFileDialogs : public QWidget
{
  Q_OBJECT
public:

  QFileDialogs( QString landmarks = "", QWidget *myParent = 0 ) : QWidget( myParent )
  {
    m_landmarks = landmarks;
  }

  QString openLandmarksFile();

  QString saveLandmarksFile();

  QString landmarksFile()
  {
    return m_landmarks;
  }

public slots:

  void openLandmarksFileSlot();

  void saveLandmarksFileSlot();

protected:

  QString m_landmarks;
};

#endif
