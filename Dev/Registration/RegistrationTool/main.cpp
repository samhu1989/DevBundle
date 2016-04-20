#include <unistd.h>
#include <iostream>
#include <QApplication>
#include "mainwindow.h"
#include <MeshViewerWidget.h>
#include <QLibraryInfo>
#include <QDebug>
#include "common.h"
#include "convertmain.h"
#include "mergemain.h"
int print_usage(int argc, char *argv[] )
{
    std::cerr << "Usage:"<<std::endl
              <<argv[0]<<" [-h] -m w " << std::endl;
    std::cerr << "Options" <<std::endl
              << "  -h\n"
              << "  Print this message\n"
              << "  -m :c Use console mode\n"
              << "      w Use ui mode.\n"
              << "  under console mode:\n"
              << "  -g : generate scene\n"
              << "  -p : convert to ply\n"
              << "  -i : input path with -g\n"
              << "     : input vertex with -p"
              << "  -n : input vertex normal\n"
              << "  -c : input vertex color\n"
              << "  -o : output file\n";

   return 0;
}


int uiMain(int argc, char *argv[])
{
    QApplication::setColorSpec( QApplication::CustomColor );
    //in case the qt is not installed
    QApplication::addLibraryPath("./qtplugins/");
    QApplication::addLibraryPath("./bin/qtplugins/");

    QApplication app(argc,argv);

    if ( !QGLFormat::hasOpenGL() ) {
      QString msg = "System has no OpenGL support!";
      QMessageBox::critical( 0, QString("OpenGL"), msg + QString(argv[1]) );
      return -1;
    }

    MainWindow mainWin;
    mainWin.resize(640,480);
    mainWin.show();

    return app.exec();
}

int cMain(int argc, char *argv[])
{
    int ch;
    opterr=0;
    while( ( ch = getopt(argc,argv,"pg") ) != -1 )
    {
        switch(ch)
        {
        case 'p':
            return convertMain(argc,argv);
        case 'g':
            return mergeMain(argc,argv);
        default:
            ;
        }
    }
}

int main(int argc, char *argv[])
{
    int ch;
    opterr=0;
    bool worked(false);
    while( ( ch = getopt(argc,argv,"hm:") ) !=-1 )
    {
        switch(ch)
        {
        case 'h':
            worked = true;
            return print_usage(argc,argv);
        case 'm':
            if(optarg[0]=='c')return cMain(argc,argv);
            else if(optarg[0]=='w') return uiMain(argc,argv);
        default:
            ;
        }
    }
    if(!worked)
    {
        return print_usage(argc,argv);
    }
    return 0;
}

