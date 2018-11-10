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

#ifndef _QLabelList_H
#define _QLabelList_H

#include "QFileDialogs.h"

#include <QFileInfo>
#include <QFile>
#include <QTextStream>
#include <QListWidget>
#include <QObject>
#include <QBrush>
#include <QColor>
#include <QPalette>

#include <cassert>
#include <iostream>
#include <map>
#include <vector>
#include <QDebug>

class QLabelList : public QListWidget
{
  Q_OBJECT

  using LandmarksMapType = std::map<QString, std::vector<double> >;
public:

  QLabelList( QWidget *myParent = nullptr ) :
    QListWidget( myParent )
  {
    m_color = 0;
  }

  void SetInputVolume( std::string filename )
  {
    m_inputVolume = QString::fromStdString( filename );
  }

  void SetInputLandmarks( std::string filename )
  {
    m_inputLandmarks = QString::fromStdString( filename );
  }

  void SetOutputLandmarks( std::string filename )
  {
    m_outputLandmarks = QString::fromStdString( filename );
  }

  LandmarksMapType GetLandmarks()
  {
    return m_landmarks;
  }

  void createListItem( const QString & label, const QString & name ); // UI for

  // creating
  // named
  // points
  // load landmarks from file
  void loadLandmarks();

  // read landmarks from label list class to the internal map
  void readLandmarks();

  // write landmarks to file
  void writeLandmarks();

public slots:

  // respond to keyboard signal from viewer
  void createListItemSlot( const QString & label );

  void createListItemAddButtonSlot(); // a wrap for add label button

  void editListItemSlot( const QString & );

  void switchListItemSlot();

  void deleteListItemSlot();

  void deleteListSlot();

  // respond to mouse signal from itself
  void cancelHighlight( QListWidgetItem * ); // haven't find a direct way of

  // disabling highlight;(

  void deleteListItemMouseSlot( QListWidgetItem * );

  void sliceChangedSlot();

  // determine which label should be displayed
  void checkVisibilitySlot();

  void checkVisibilitySlot( double *tag ); // a wrap for slider bar

  void checkVisibilitySlot( QListWidgetItem * ); // a wrap for double click on

  // list item

  // respond to callback, help to find the initial position due to wheeling
  void ackWheelChanged();

  // save landmarks to a file
  void saveLandmarks();

  void saveAsLandmarks();

signals:

  void sliceChangedList(); // signal to sliceChangeSlot indicating the update of

  // slice viewer

  void sendLabelPosition( double *pos );

  // determine which label should be displayed
  // *[0] = sagittal, *[1] = coronal, *[2] = axial
  // table = 1 means the label is visible in certain viewer
  void visibilityTable( int *table );

  // send to viewer, help to find the initial position due to wheeling
  void sendLabelPositions( double *pos );

protected:

  // color seed
  int m_color;

  // label postion signal received
  QString m_label;

  LandmarksMapType m_landmarks;

  QString m_inputVolume;
  QString m_inputLandmarks;
  QString m_outputLandmarks;
};

#endif
