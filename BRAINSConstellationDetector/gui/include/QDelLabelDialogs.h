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

class QDelLabelDialogs : public QDialog
{
  Q_OBJECT
public:

  QDelLabelDialogs( QString text = "", QWidget *myParent = 0 );
public slots:

  void exec( QListWidgetItem * );

  void accept2(); // a wrap for double click situation

signals:

  void exec2(); // a wrap for double click situation

  void accepted( QListWidgetItem * );

protected:

  QGridLayout *    m_layout;
  QLabel           m_label;
  QPushButton      m_accept;
  QPushButton      m_cancel;
  QListWidgetItem *m_listItem;
};

#endif
