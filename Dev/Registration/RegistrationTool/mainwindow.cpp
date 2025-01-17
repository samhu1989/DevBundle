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
    opt += OpenMesh::IO::Options::Binary;
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
    connect(ui->actionJRMPC,SIGNAL(triggered(bool)),this,SLOT(start_registration()));
    connect(ui->actionPMSDP,SIGNAL(triggered(bool)),this,SLOT(start_registration()));

    connect(&timer,SIGNAL(timeout()),w,SLOT(updateGL()));
    timer.setSingleShot(false);
    //force the repaint in gl every 100ms
    timer.start(50);
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
        Config::Ptr config_ptr;
        if(!thread->init(w->list(),config_ptr))
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
    if(s==ui->actionPMSDP)
    {
        PMSDP_Thread* thread = new PMSDP_Thread();
        Config::Ptr config_ptr;
        if(!thread->init(w->list(),config_ptr))
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
    if(s==ui->actionJRMPC)
    {
        JRMPC_Thread* thread = new JRMPC_Thread();
        Config::Ptr config_ptr;
        if(!thread->init(w->list(),config_ptr))
        {
            QString msg = "Fail to Initialize the Registration:\n '";
            msg += QString::fromStdString(thread->errorString());
            QMessageBox::critical( NULL, windowTitle(), msg);
            return;
        }else{
            w->reset_center();
            w->show_back();
            connect(thread,SIGNAL(finished()),this,SLOT(finish_registration()));
        }
        alg_thread = thread;
    }
    if(alg_thread){
        time.restart();
        alg_thread->start(QThread::HighestPriority);
    }
}

void MainWindow::finish_registration(void)
{
    std::cerr<<time.elapsed()<<std::endl;
    QString msg ;
    msg = msg.sprintf("Current Algorithm is Finished(%d ms)",time.elapsed()) ;
    if(alg_thread)
    {
        while(alg_thread->isRunning())
        {
            alg_thread->terminate();
            QApplication::processEvents();
        }
        alg_thread->deleteLater();
        alg_thread = NULL;
    }

    QMessageBox::information( this, windowTitle(), msg);
    ui->statusbar->showMessage(msg);
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
            alg_thread->quit();
        }
    }
    QMainWindow::keyPressEvent(e);
}

MainWindow::~MainWindow()
{
    delete ui;
}
