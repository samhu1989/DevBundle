#include "mainwindow.h"
#include "ui_mainwindow.h"
#include "MeshPairViewerWidget.h"
#include "computenormalthread.h"
#include "computeoctreethread.h"
#include "computesupervoxelthread.h"
#include "downsamplethread.h"
#include "QMessageBox"
#include "extractplanethread.h"
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
    MeshBundle<DefaultMesh>::Ptr output(new MeshBundle<DefaultMesh>);
    ComputeOctreeThread th(v->first_ptr()->mesh_,output->mesh_);
    th.start(QThread::HighestPriority);
    while(!th.wait(30))
    {
        QApplication::processEvents();
    }
    v->second_ptr() = output;
    QString msg = "Done Compute Octree:\n '";
    QMessageBox::information( this, windowTitle(), msg );
}

void MainWindow::computeSuperVoxel()
{
    MeshPairViewerWidget* v = qobject_cast<MeshPairViewerWidget*>(centralWidget());
    if(!v)return;
    ComputeSupervoxelThread th(v->first_ptr(),v->second_ptr());
    th.start(QThread::HighestPriority);
    while(!th.wait(30))
    {
        QApplication::processEvents();
    }
    QString msg = "Done Compute SuperVoxel:\n '";
    QMessageBox::information( this, windowTitle(), msg );
}

void MainWindow::computeDownSample()
{
    MeshPairViewerWidget* v = qobject_cast<MeshPairViewerWidget*>(centralWidget());
    if(!v)return;
    DownSampleThread th(v->first_ptr());
    th.start(QThread::HighestPriority);
    while(!th.wait(30))
    {
        QApplication::processEvents();
    }
    QString msg = "Done Down Sample:\n '";
    QMessageBox::information( this, windowTitle(), msg );
}

void MainWindow::computeExtractPlane()
{
    MeshPairViewerWidget* v = qobject_cast<MeshPairViewerWidget*>(centralWidget());
    if(!v)return;
    ExtractPlaneThread th(v->first_ptr());
    th.start(QThread::HighestPriority);
    while(!th.wait(30))
    {
        QApplication::processEvents();
    }
    QString msg = "Done Extract Plane:\n '";
    QMessageBox::information( this, windowTitle(), msg );
}

MainWindow::~MainWindow()
{
    delete ui;
}
