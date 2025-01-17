/* ========================================================================= *
 *                                                                           *
 *                               OpenMesh                                    *
 *           Copyright (c) 2001-2015, RWTH-Aachen University                 *
 *           Department of Computer Graphics and Multimedia                  *
 *                          All rights reserved.                             *
 *                            www.openmesh.org                               *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * This file is part of OpenMesh.                                            *
 *---------------------------------------------------------------------------*
 *                                                                           *
 * Redistribution and use in source and binary forms, with or without        *
 * modification, are permitted provided that the following conditions        *
 * are met:                                                                  *
 *                                                                           *
 * 1. Redistributions of source code must retain the above copyright notice, *
 *    this list of conditions and the following disclaimer.                  *
 *                                                                           *
 * 2. Redistributions in binary form must reproduce the above copyright      *
 *    notice, this list of conditions and the following disclaimer in the    *
 *    documentation and/or other materials provided with the distribution.   *
 *                                                                           *
 * 3. Neither the name of the copyright holder nor the names of its          *
 *    contributors may be used to endorse or promote products derived from   *
 *    this software without specific prior written permission.               *
 *                                                                           *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS       *
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED *
 * TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A           *
 * PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER *
 * OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,  *
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,       *
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR        *
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF    *
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING      *
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS        *
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.              *
 *                                                                           *
 * ========================================================================= */

/*===========================================================================*\
 *                                                                           *
 *   $Revision: 1258 $                                                         *
 *   $Date: 2015-04-28 15:07:46 +0200 (Di, 28 Apr 2015) $                   *
 *                                                                           *
\*===========================================================================*/

#ifdef _MSC_VER
#  pragma warning(disable: 4267 4311)
#endif

#include <iostream>
#include <fstream>
#include <QApplication>
#include <QMessageBox>
#include <QMainWindow>
#include <QMenuBar>
#include <QFileDialog>
#ifdef ARCH_DARWIN
#include <glut.h>
#else
#include <GL/glut.h>
#endif
#include "MeshColor.h"
#include "MeshPairViewerWidget.h"
#include <iomanip>
#include "mainwindow.h"
#include <QTimer>

void create_menu(MainWindow &w);
void usage_and_exit(int xcode);

int main(int argc, char **argv)
{
  // OpenGL check
  QApplication::setColorSpec( QApplication::CustomColor );
  QApplication::addLibraryPath("./qtplugins/");
  QApplication::addLibraryPath("./bin/qtplugins/");
  QApplication app(argc,argv);
#if !defined(__APPLE__)
  glutInit(&argc,argv);
#endif

  if ( !QGLFormat::hasOpenGL() ) {
    QString msg = "System has no OpenGL support!";
    QMessageBox::critical( 0, QString("OpenGL"), msg + QString(argv[1]) );
    return -1;
  }

  int c;
  OpenMesh::IO::Options opt;

  while ( (c=getopt(argc,argv,"hbs"))!=-1 )
  {
     switch(c)
     {
       case 'b': opt += OpenMesh::IO::Options::Binary; break;
       case 'h':
          usage_and_exit(0);
       case 's': opt += OpenMesh::IO::Options::Swap; break;
       default:
          usage_and_exit(1);
     }
  }

  // enable most options for now
  opt += OpenMesh::IO::Options::VertexColor;
  opt += OpenMesh::IO::Options::VertexNormal;
//  opt += OpenMesh::IO::Options::VertexTexCoord;
//  opt += OpenMesh::IO::Options::FaceColor;
//  opt += OpenMesh::IO::Options::FaceNormal;
//  opt += OpenMesh::IO::Options::FaceTexCoord;

  // create widget
  MainWindow mainWin;
  MeshPairViewerWidget w(&mainWin);
  QTimer gl_timer;
  gl_timer.setSingleShot(false);
  QObject::connect(&gl_timer,SIGNAL(timeout()),&w,SLOT(updateGL()));
  w.setOptions(opt);
  mainWin.setCentralWidget(&w);

  create_menu(mainWin);

  // static mesh, hence use strips
  w.enable_strips();

  mainWin.resize(640, 480);
  mainWin.show();

  // load scene if specified on the command line
//  if ( optind < argc )
//  {
//    w.open_mesh_gui(argv[optind]);
//  }

//  if ( ++optind < argc )
//  {
//      w.open_texture_gui(argv[optind]);
//  }

  gl_timer.start(100);
  return app.exec();
}

void create_menu(MainWindow &w)
{
    using namespace Qt;
    QMenu *fileMenu = w.menuBar()->addMenu(w.tr("&File"));

    QAction* openAct = new QAction(w.tr("&Open mesh..."), &w);
    openAct->setShortcut(w.tr("Ctrl+O"));
    openAct->setStatusTip(w.tr("Open a mesh file"));
    QObject::connect(openAct, SIGNAL(triggered()), w.centralWidget(), SLOT(query_open_source_file()));
    fileMenu->addAction(openAct);

    QAction* saveAct = new QAction(w.tr("&Save mesh..."), &w);
    saveAct->setShortcut(w.tr("Ctrl+S"));
    saveAct->setStatusTip(w.tr("Save a mesh file"));
    QObject::connect(saveAct, SIGNAL(triggered()), w.centralWidget(), SLOT(query_save_mesh_file()));
    fileMenu->addAction(saveAct);

    QMenu *editMenu = w.menuBar()->addMenu(w.tr("&Edit"));

    QAction* normalAct = new QAction(w.tr("&Compute Vertex Normal"), &w);
    normalAct->setStatusTip(w.tr("Compute Vertex Normal"));
    QObject::connect(normalAct, SIGNAL(triggered()), &w, SLOT(computeVertexNormal()));
    editMenu->addAction(normalAct);

    QAction* octreeAct = new QAction(w.tr("&Compute Octree"), &w);
    octreeAct->setStatusTip(w.tr("Compute Octree"));
    QObject::connect(octreeAct, SIGNAL(triggered()), &w, SLOT(computeOctree()));
    editMenu->addAction(octreeAct);

    QAction* svAct = new QAction(w.tr("&Compute Supervoxel"), &w);
    svAct->setStatusTip(w.tr("Compute Supervoxel"));
    QObject::connect(svAct, SIGNAL(triggered()), &w, SLOT(computeSuperVoxel()));
    editMenu->addAction(svAct);

    QAction* dsAct = new QAction(w.tr("&Down Sample"), &w);
    dsAct->setStatusTip(w.tr("Down Sample"));
    QObject::connect(dsAct, SIGNAL(triggered()), &w, SLOT(computeDownSample()));
    editMenu->addAction(dsAct);

    QAction* epAct = new QAction(w.tr("&Extract Plane"), &w);
    epAct->setStatusTip(w.tr("Extract Plane"));
    QObject::connect(epAct, SIGNAL(triggered()), &w, SLOT(computeExtractPlane()));
    editMenu->addAction(epAct);

    QMenu *viewMenu = w.menuBar()->addMenu(w.tr("&View"));

    QAction* customcolorAct = new QAction(w.tr("&Show Custom Color"), &w);
    customcolorAct->setCheckable(true);
    customcolorAct->setStatusTip(w.tr("Show Custom Color"));
    QObject::connect(customcolorAct, SIGNAL(toggled(bool)), w.centralWidget(), SLOT(use_custom_color(bool)));
    viewMenu->addAction(customcolorAct);

}

void usage_and_exit(int xcode)
{
   std::cout << "Usage: meshviewer [-s] [mesh] [texture]\n" << std::endl;
   std::cout << "Options:\n"
	     << "  -b\n"
	     << "    Assume input to be binary.\n\n"
             << "  -s\n"
             << "    Reverse byte order, when reading binary files.\n"
             << std::endl;
   exit(xcode);
}
