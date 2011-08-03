/*
 * Author: Wei Lu
 * at Psychiatry Imaging Lab,
 * University of Iowa Health Care 2010
 */

#ifndef __QHelpDialog_H
#define __QHelpDialog_H

#include <QDialog>
#include <QGridLayout>
#include <QLabel>
#include <QObject>

class QHelpDialog : public QDialog
{
  Q_OBJECT
public:

  QHelpDialog( QWidget *myParent = 0 );
protected:

  QGridLayout *m_layout;
  QLabel       m_label;
};

#endif
