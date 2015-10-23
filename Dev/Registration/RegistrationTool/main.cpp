#include <unistd.h>
#include <iostream>
#include <QApplication>
#include "mainwindow.h"
#include <MeshViewerWidget.hh>

int print_usage(int argc, char *argv[] )
{
    std::cerr << "Usage:"<<std::endl
              <<argv[0]<<" [-h] -m w " << std::endl;
    std::cerr << "Options" <<std::endl
              << "  -h\n"
              << "  Print this message\n"
              << "  -m\n"
              << "  c : Use console mode\n"
              << "  else : Use ui mode.\n";
   return 0;
}


int uiMain(int argc, char *argv[])
{
    QApplication::setColorSpec( QApplication::CustomColor );
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
    std::cerr<<"cMain"<<std::endl;
    return 0;
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
            else return uiMain(argc,argv);
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

