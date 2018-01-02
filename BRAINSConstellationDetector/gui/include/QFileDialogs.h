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

#ifndef _QFileDialogs_H
#define _QFileDialogs_H

#include <QFileDialog>
#include <QObject>

#include "itkMacro.h" //Needed for nullptr

class QFileDialogs : public QWidget
{
  Q_OBJECT
public:

  QFileDialogs( QString landmarks = "", QWidget *myParent = nullptr ) : QWidget( myParent )
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
