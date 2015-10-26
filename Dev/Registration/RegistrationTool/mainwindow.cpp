#include "mainwindow.h"
#include "ui_mainwindow.h"
#include "MeshListViewerWidget.h"
MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);
    MeshListViewerWidget* w = new MeshListViewerWidget(this);
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
    connect(ui->actionOpen,SIGNAL(triggered(bool)),w,SLOT(query_open_file()));
    connect(ui->actionSave,SIGNAL(triggered(bool)),w,SLOT(query_save_file()));
}

MainWindow::~MainWindow()
{
    delete ui;
}
