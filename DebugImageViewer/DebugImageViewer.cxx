#include <QApplication>

#include "QDebugImageViewerWindow.h"

int main(int argc, char *argv[])
{
  QApplication            app(argc, argv);
  QDebugImageViewerWindow window;
  QObject::connect(&app, SIGNAL(aboutToQuit() ),
                   &window, SLOT(exiting() ) );

  window.show();
  exit(app.exec() );
}
