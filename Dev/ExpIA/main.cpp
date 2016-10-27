#include "mainwindow.h"
#include <QApplication>
#include "common.h"
int main(int argc, char *argv[])
{
    QApplication::setColorSpec( QApplication::CustomColor );
    QApplication::addLibraryPath("./qtplugins/");
    QApplication::addLibraryPath("./bin/qtplugins/");
    init_resouce();
    QApplication a(argc, argv);
    MainWindow w;
    w.show();
    return a.exec();
}
