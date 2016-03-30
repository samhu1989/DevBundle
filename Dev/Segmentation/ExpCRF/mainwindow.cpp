#include "mainwindow.h"
#include "ui_mainwindow.h"
#include <QMessageBox>
#include "crf2d.h"
MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    edit_thread_(NULL),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);
    connect(ui->actionCRF,SIGNAL(triggered(bool)),this,SLOT(start_editing()));
    connect(ui->actionCRF3D,SIGNAL(triggered(bool)),this,SLOT(start_editing()));
}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::start_editing()
{
    if(edit_thread_)
    {
        QString msg = "Please Wait Till the End of Last Algorithm\n";
        QMessageBox::critical(this, windowTitle(), msg);
        return;
    }
    QAction* edit = qobject_cast<QAction*>(sender());
    if(edit==ui->actionCRF)
    {
        edit_thread_ = new QThread();
        CRF2D* worker = new CRF2D();
        worker->moveToThread(edit_thread_);
        connect(edit_thread_,SIGNAL(started()),worker,SLOT(process()));
    }
    connect(edit_thread_,SIGNAL(finished()),this,SLOT(finish_editing()));
    edit_thread_->start(QThread::HighPriority);
}

void MainWindow::finish_editing()
{
    if(edit_thread_)
    {
        edit_thread_->deleteLater();
        edit_thread_ = NULL;
    }
}
