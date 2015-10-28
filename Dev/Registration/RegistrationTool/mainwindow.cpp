#include "mainwindow.h"
#include "ui_mainwindow.h"
#include "MeshListViewerWidget.h"
#include "registrationcore.h"
MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow),
    alg_thread(NULL)
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
    connect(ui->actionCPDRigid3D,SIGNAL(triggered(bool)),this,SLOT(start_registration()));
    connect(&timer,SIGNAL(timeout()),w,SLOT(updateGL()));
    timer.setSingleShot(false);
    //force the repaint in gl every 100ms
    timer.start(100);
}

void MainWindow::start_registration(void)
{
    typedef std::vector<MeshBundle<DefaultMesh>::Ptr> MeshList;
    if(alg_thread)
    {
        QString msg = "Please Wait Till the End of Last Algorithm:\n '";
        QMessageBox::critical( NULL, windowTitle(), msg);
        return;
    }
    QAction* s = qobject_cast<QAction*>(sender());
    MeshListViewerWidget* w;
    if(!s)return;
    else{
        w = qobject_cast<MeshListViewerWidget*>(centralWidget());
        if(!w){
            QString msg = "No Visualization Port:\n '";
            QMessageBox::critical( NULL, windowTitle(), msg);
            return;
        }
    }
    if(s==ui->actionCPDRigid3D)
    {
        CPDR3D_DM_R_Thread* thread = new CPDR3D_DM_R_Thread();
        if(!thread->init(w->list()))
        {
            QString msg = "Fail to Initialize the Registration:\n '";
            msg += QString::fromStdString(thread->errorString());
            QMessageBox::critical( NULL, windowTitle(), msg);
            return;
        }else{
            connect(thread,SIGNAL(finished()),this,SLOT(finish_registration()));
        }
        alg_thread = thread;
    }
    if(alg_thread)alg_thread->start(QThread::NormalPriority);
}

void MainWindow::finish_registration(void)
{
    QString msg = "Current Algorithm is finished";
    QMessageBox::information( this, windowTitle(), msg);
    if(alg_thread)
    {
        if(!alg_thread->isRunning()){
            alg_thread->deleteLater();
            alg_thread = NULL;
        }
    }
}

void MainWindow::keyPressEvent(QKeyEvent*e)
{
    if(e->key()==Qt::Key_K)
    {
        QString msg = "Do you really want to kill current registration?\n '";
        if(QMessageBox::Yes == QMessageBox::information(
                    NULL,
                    windowTitle(),
                    msg,
                    QMessageBox::Yes|QMessageBox::No,
                    QMessageBox::No)
                )
        {
            alg_thread->terminate();
        }
    }
    QMainWindow::keyPressEvent(e);
}

MainWindow::~MainWindow()
{
    delete ui;
}
