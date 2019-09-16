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

#ifndef __QDelLabelDialogs_H
#define __QDelLabelDialogs_H

#include <QDialog>
#include <QGridLayout>
#include <QLabel>
#include <QPushButton>
#include <QListWidget>
#include <QObject>

#include "itkMacro.h" //Needed for nullptr

class QDelLabelDialogs : public QDialog
{
  Q_OBJECT
public:
  QDelLabelDialogs(QString text = "", QWidget * myParent = nullptr);
public slots:

  void
  exec(QListWidgetItem *);

  void
  accept2(); // a wrap for double click situation

signals:

  void
  exec2(); // a wrap for double click situation

  void
  accepted(QListWidgetItem *);

protected:
  QGridLayout *     m_layout;
  QLabel            m_label;
  QPushButton       m_accept;
  QPushButton       m_cancel;
  QListWidgetItem * m_listItem;
};

#endif
