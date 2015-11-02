#include "mainwindow.h"
#include "ui_mainwindow.h"
#include "MeshViewerWidget.h"
#include "computenormalthread.h"
#include "QMessageBox"
MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);
}

void MainWindow::computeVertexNormal()
{
    MeshViewerWidget* v = qobject_cast<MeshViewerWidget*>(centralWidget());
    if(!v)return;
    ComputeNormalThread th(v->mesh());
    th.start(QThread::HighestPriority);
    while(!th.wait(30))
    {
        QApplication::processEvents();
    }
    QString msg = "Done Compute Normal:\n '";
    QMessageBox::information( this, windowTitle(), msg );
}

MainWindow::~MainWindow()
{
    delete ui;
}
