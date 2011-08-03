/*
 * Author: Wei Lu
 * at Psychiatry Imaging Lab,
 * University of Iowa Health Care 2010
 */

#include "QFileDialogs.h"

QString QFileDialogs::openLandmarksFile()
{
  m_landmarks =
    QFileDialog::getOpenFileName( this,
                                  tr("Select landmarks file"),
                                  QDir::currentPath(),
                                  tr("Model file ( *.fcsv )") );
  return m_landmarks;
}

void QFileDialogs::openLandmarksFileSlot()
{
  m_landmarks =
    QFileDialog::getOpenFileName( this,
                                  tr("Select landmarks file"),
                                  QDir::currentPath(),
                                  tr("Model file ( *.fcsv )") );
}

QString QFileDialogs::saveLandmarksFile()
{
  m_landmarks =
    QFileDialog::getSaveFileName( this,
                                  tr("Save landmarks file"),
                                  QDir::currentPath(),
                                  tr("Model file ( *.fcsv )") );
  return m_landmarks;
}

void QFileDialogs::saveLandmarksFileSlot()
{
  m_landmarks =
    QFileDialog::getSaveFileName( this,
                                  tr("Save landmarks file"),
                                  QDir::currentPath(),
                                  tr("Model file ( *.fcsv )") );
}
