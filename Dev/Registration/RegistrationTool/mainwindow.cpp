#include "mainwindow.h"
#include "ui_mainwindow.h"
#include "MeshViewerWidget.hh"
MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);
    MeshViewerWidget* w = new MeshViewerWidget(this);
    OpenMesh::IO::Options opt;
    opt += OpenMesh::IO::Options::VertexColor;
    opt += OpenMesh::IO::Options::VertexNormal;
    opt += OpenMesh::IO::Options::VertexTexCoord;
    opt += OpenMesh::IO::Options::FaceColor;
    opt += OpenMesh::IO::Options::FaceNormal;
    opt += OpenMesh::IO::Options::FaceTexCoord;
    w->setOptions(opt);
    setCentralWidget(w);
    connect(this,SIGNAL(destroyed(QObject*)),w,SLOT(deleteLater()));
    connect(ui->actionOpen_Source,SIGNAL(triggered(bool)),w,SLOT(query_open_mesh_file()));
}

MainWindow::~MainWindow()
{
    delete ui;
}
