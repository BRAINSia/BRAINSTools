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
#include <itkSpatialOrientationAdapter.h>
#include "QDebugImageViewerWindow.h"
#include "QImageDisplay.h"
#include <string>
#include <QApplication>
#include <QString>
#include <QHBoxLayout>
#include <QTcpServer>
#include <QTcpSocket>

QDebugImageViewerWindow::QDebugImageViewerWindow(QWidget * parent, Qt::WindowFlags flags)
  : QMainWindow(parent, flags)
{
  this->m_ViewCount = 1;
  std::string cmdLineImageName;
  for (int i = 0; i < qApp->arguments().size(); i++)
  {
    QString current(qApp->arguments().at(i));
    if (current == "--numviews")
    {
      ++i;
      if (i >= qApp->arguments().size())
      {
        break;
      }
      QString currentArg(qApp->arguments().at(i));
      this->m_ViewCount = currentArg.toInt();
    }
    else if (current == "--input")
    {
      ++i;
      if (i >= qApp->arguments().size())
      {
        break;
      }
      QString currentArg(qApp->arguments().at(i));
      cmdLineImageName = currentArg.toStdString();
    }
  }
  QWidget * main = new QWidget(this);
  main->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);
  this->setCentralWidget(main);

  QHBoxLayout * mainLayout = new QHBoxLayout(main);
  main->setLayout(mainLayout);
  for (int i = 0; i < this->m_ViewCount; i++)
  {
    QImageDisplay * current = new QImageDisplay(this);
    mainLayout->addWidget(current);
    this->m_ImageDisplayList.push_back(current);
    if (cmdLineImageName != "")
    {
      current->SetImage(cmdLineImageName);
    }
    else
    {
      current->SetBlankImage();
    }
  }
  // set up server
  this->m_Server = new QTcpServer(this);
  connect(this->m_Server, SIGNAL(newConnection()), this, SLOT(newConnection()));
  if (!this->m_Server->listen(QHostAddress::LocalHost, 19345))
  {
    std::cerr << "Can't start server on port 19345" << std::endl;
  }
}

void
QDebugImageViewerWindow::newConnection()
{
  this->m_Socket = this->m_Server->nextPendingConnection();
  std::cerr << "buffer size " << this->m_Socket->readBufferSize() << std::endl;
  this->SetupSocketConnections();
}

void
QDebugImageViewerWindow::SetupSocketConnections()
{
  connect(this->m_Socket, SIGNAL(readyRead()), this, SLOT(readImage()));

  connect(this->m_Socket,
          SIGNAL(stateChanged(QAbstractSocket::SocketState)),
          this,
          SLOT(stateChanged(QAbstractSocket::SocketState)));
}

//
// this function sort of apes the fread signature
qint64
QDebugImageViewerWindow::SocketRead(void * buf, qint64 objectSize, qint64 objectCount)
{
  qint64 toRead = objectSize * objectCount;

  qint64 readCount(0);
  qint64 leftToRead(toRead);

  char * readPtr(reinterpret_cast<char *>(buf));

  while (leftToRead > 0)
  {
    this->m_Socket->waitForReadyRead();
    qint64 readAmt = this->m_Socket->read(readPtr, leftToRead);

    switch (readAmt)
    {
      case -1:
        return -1;
      case 0:
        return readCount;
      default:
        readPtr += readAmt;
        readCount += readAmt;
        leftToRead -= readAmt;
    }
  }

  return readCount / objectSize;
}

void
QDebugImageViewerWindow::stateChanged(QAbstractSocket::SocketState state)
{
  switch (state)
  {
    //   6       The socket is about to close (data may
    //           still be waiting to be written).
    case QAbstractSocket::ClosingState:
    //   0       The socket is not connected.
    case QAbstractSocket::UnconnectedState:
      break;
    //   1       The socket is performing a host name lookup.
    case QAbstractSocket::HostLookupState:
      break;
    //   2       The socket has started establishing a connection.
    case QAbstractSocket::ConnectingState:
      break;
    //   3       A connection is established.
    case QAbstractSocket::ConnectedState:
      break;
    //   4       The socket is bound to an address and port (for servers).
    case QAbstractSocket::BoundState:
      break;
    //   5       For internal use only.
    case QAbstractSocket::ListeningState:
      break;
  }
}

void
QDebugImageViewerWindow::readImage()
{
  QImageDisplay::ImageType::SizeType imageSize;

  for (unsigned i = 0; i < 3; i++)
  {
    if (this->SocketRead(&imageSize[i], sizeof(imageSize[0]), 1) != 1)
    {
      std::cerr << "Error reading socket" << std::endl;
      exit(1);
      return;
    }
  }

  QImageDisplay::ImageType::SpacingType imageSpacing;
  QImageDisplay::ImageType::IndexType   imageIndex;
  for (unsigned i = 0; i < 3; i++)
  {
    imageIndex[i] = 0; // sneak initializing index into this loop
    if (this->SocketRead(&imageSpacing[i], sizeof(imageSpacing[0]), 1) != 1)
    {
      std::cerr << "Error reading socket" << std::endl;
      exit(1);
    }
  }
  QImageDisplay::ImageType::RegionType imageRegion;
  imageRegion.SetSize(imageSize);
  imageRegion.SetIndex(imageIndex);

  QImageDisplay::ImageType::Pointer xferImage = QImageDisplay::ImageType::New();
  xferImage->SetSpacing(imageSpacing);
  xferImage->SetRegions(imageRegion);
  xferImage->Allocate();

  // set orientation
  itk::SpatialOrientation::ValidCoordinateOrientationFlags orientation;
  if (this->SocketRead(&orientation, sizeof(orientation), 1) != 1)
  {
    std::cerr << "Error reading socket" << std::endl;
    exit(1);
  }
  xferImage->SetDirection(itk::SpatialOrientationAdapter().ToDirectionCosines(orientation));

  QImageDisplay::ImageType::PointType origin;
  for (unsigned i = 0; i < 3; i++)
  {
    double oVal;
    if (this->SocketRead(&oVal, sizeof(double), 1) != 1)
    {
      std::cerr << "Error reading socket" << std::endl;
      exit(1);
    }
    origin[i] = oVal;
  }
  xferImage->SetOrigin(origin);

  unsigned int viewIndex;
  if (this->SocketRead(&viewIndex, sizeof(viewIndex), 1) != 1)
  {
    std::cerr << "Error reading socket" << std::endl;
    exit(1);
  }
  unsigned char * pixelData = xferImage->GetBufferPointer();

  qint64 bufferSize(imageSize[0] * imageSize[1] * imageSize[2]);
  if (this->SocketRead(pixelData, sizeof(QImageDisplay::ImageType::PixelType), bufferSize) != bufferSize)
  {
    std::cerr << "Error reading socket" << std::endl;
    exit(1);
  }
  if (viewIndex < this->m_ViewCount)
  {
    this->m_ImageDisplayList[viewIndex]->SetImage(xferImage);
  }
}

void
QDebugImageViewerWindow::exiting()
{}
