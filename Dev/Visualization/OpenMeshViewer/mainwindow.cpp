#include "mainwindow.h"
#include "ui_mainwindow.h"
#include "MeshPairViewerWidget.h"
#include "computenormalthread.h"
#include "computeoctreethread.h"
#include "QMessageBox"
MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);
}

void MainWindow::computeVertexNormal()
{
    MeshPairViewerWidget* v = qobject_cast<MeshPairViewerWidget*>(centralWidget());
    if(!v)return;
    ComputeNormalThread th(v->first_ptr()->mesh_);
    th.start(QThread::HighestPriority);
    while(!th.wait(30))
    {
        QApplication::processEvents();
    }
    QString msg = "Done Compute Normal:\n '";
    QMessageBox::information( this, windowTitle(), msg );
}

void MainWindow::computeOctree()
{
    MeshPairViewerWidget* v = qobject_cast<MeshPairViewerWidget*>(centralWidget());
    if(!v)return;
    ComputeOctreeThread th(v->first_ptr()->mesh_,v->second_ptr()->mesh_);
    th.start(QThread::HighestPriority);
    while(!th.wait(30))
    {
        QApplication::processEvents();
    }
    QString msg = "Done Compute Octree:\n '";
    QMessageBox::information( this, windowTitle(), msg );
}

MainWindow::~MainWindow()
{
    delete ui;
}
