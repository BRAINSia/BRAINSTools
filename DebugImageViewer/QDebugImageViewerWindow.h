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
#ifndef QDebugImageViewerWindow_h
#define QDebugImageViewerWindow_h
#include <QMainWindow>
#include <QAbstractSocket> >
#include <vector>

class QImageDisplay;
class QTcpServer;
class QTcpSocket;

class QDebugImageViewerWindow : public QMainWindow
{
  Q_OBJECT
public slots:
  void
  exiting();

  void
  newConnection();

  void
  readImage();

  void
  stateChanged(QAbstractSocket::SocketState state);

public:
  QDebugImageViewerWindow(QWidget * parent = 0, Qt::WindowFlags flags = 0);

private:
  qint64
  SocketRead(void * buf, qint64 objectSize, qint64 objectCount);

  void
  SetupSocketConnections();

  using ImageDisplayListType = std::vector<QImageDisplay *>;
  ImageDisplayListType m_ImageDisplayList;
  int                  m_ViewCount;
  QTcpServer *         m_Server;
  QTcpSocket *         m_Socket;
};

#endif // QDebugImageViewerWindow_h
