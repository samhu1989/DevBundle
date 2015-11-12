#include "mainwindow.h"
#include <QApplication>

int main(int argc, char *argv[])
{
    QApplication::setColorSpec( QApplication::CustomColor );
    QApplication::addLibraryPath("./qtplugins/");
    QApplication::addLibraryPath("./bin/qtplugins/");
    QApplication a(argc, argv);
    MainWindow w;
    w.show();
    return a.exec();
}
