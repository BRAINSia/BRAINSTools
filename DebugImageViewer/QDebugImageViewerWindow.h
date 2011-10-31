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
  void exiting();

  void newConnection();

  void readImage();

  void stateChanged(QAbstractSocket::SocketState state);

public:
  QDebugImageViewerWindow(QWidget *parent = 0,
                          Qt::WindowFlags flags = 0);
private:
  qint64 SocketRead(void *buf, qint64 objectSize, qint64 objectCount);

  void SetupSocketConnections();

  typedef std::vector<QImageDisplay *> ImageDisplayListType;
  ImageDisplayListType m_ImageDisplayList;
  int                  m_ViewCount;
  QTcpServer*          m_Server;
  QTcpSocket*          m_Socket;
};

#endif // QDebugImageViewerWindow_h
