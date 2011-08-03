/*
 * Author: Wei Lu
 * at Psychiatry Imaging Lab,
 * University of Iowa Health Care 2010
 */

#include "QDelLabelDialogs.h"

QDelLabelDialogs::QDelLabelDialogs(QString text, QWidget *myParent) : QDialog(myParent)
{
  m_label.setText(text);
  m_accept.setText("Accept");
  m_cancel.setText("Cancel");

  m_layout = new QGridLayout(this);
  m_layout->addWidget(&m_label, 1, 1, Qt::AlignHCenter);
  m_layout->addWidget(&m_accept, 2, 1, Qt::AlignRight);
  m_layout->addWidget(&m_cancel, 2, 2, Qt::AlignRight);

  // for buttons and keyboard input
  QObject::connect( &m_accept, SIGNAL( clicked() ),
                    this, SLOT( accept() ) );
  QObject::connect( &m_cancel, SIGNAL( clicked() ),
                    this, SLOT( reject() ) );

  // for mouse double click
  QObject::connect( this, SIGNAL( exec2() ),
                    this, SLOT( exec() ) );
  QObject::connect( &m_accept, SIGNAL( clicked() ),
                    this, SLOT( accept2() ) );

  m_listItem = NULL;
}

void QDelLabelDialogs::exec(QListWidgetItem *listItem)
{
  m_listItem = listItem;
  emit exec2();
}

void QDelLabelDialogs::accept2()
{
  emit accepted(m_listItem);
}
