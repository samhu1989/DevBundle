#include "mainwindow.h"
#include "ui_mainwindow.h"
#include <QMessageBox>
#include "crf2d.h"
#include <QFileDialog>
#include <QLabel>
MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    edit_thread_(NULL),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);
    connect(ui->actionCRF,SIGNAL(triggered(bool)),this,SLOT(start_editing()));
    connect(ui->actionCRF3D,SIGNAL(triggered(bool)),this,SLOT(start_editing()));
    connect(ui->actionLoad_Image,SIGNAL(triggered(bool)),this,SLOT(load_img()));
}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::load_img(void)
{
    QString fileName = QFileDialog::getOpenFileName(this,
        tr("Load Input Image"),
        tr("../Dev_Data"),
        tr(";;"
        "Images Files (*.ppm,*.png,*jpg)"));
    if (fileName.isEmpty())return;

    if(!input_img_.load(fileName))
    {
        QString msg = "Failed to load "+fileName+"\n";
        QMessageBox::critical(this, windowTitle(), msg);
        return;
    }
    QLabel* imgview = new QLabel();
    imgview->setAttribute(Qt::WA_DeleteOnClose,true);
    imgview->setPixmap(QPixmap::fromImage(input_img_));
    imgview->setAlignment(Qt::AlignCenter|Qt::AlignHCenter);
    ui->mdiArea->addSubWindow(imgview);
    imgview->show();
}

void MainWindow::load_annotation(void)
{
    QString fileName = QFileDialog::getOpenFileName(this,
        tr("Load Input Image"),
        tr("./data"),
        tr(";;"
        "All Files (*.arma,*.ppm)"));
    if (!fileName.isEmpty())
    {
        if(input_img_.load(fileName))return;
        QString msg = "Failed to load "+fileName+"\n";
        QMessageBox::critical(this, windowTitle(), msg);
    }
}

void MainWindow::start_editing()
{
    if(edit_thread_)
    {
        QString msg = "Please Wait Till the End of Last Algorithm\n";
        QMessageBox::critical(this, windowTitle(), msg);
        return;
    }
    if(input_img_.isNull())
    {
        QString msg = "Please Load Input First\n";
        QMessageBox::critical(this, windowTitle(), msg);
        return;
    }
    if(annotation_.empty())
    {
        QString msg = "Please Load Annotation First\n";
        QMessageBox::critical(this, windowTitle(), msg);
        return;
    }
    QAction* edit = qobject_cast<QAction*>(sender());
    if(edit==ui->actionCRF)
    {
        edit_thread_ = new QThread();
        CRF2D* worker = new CRF2D(input_img_,annotation_);
        connect(worker,SIGNAL(message(QString,int)),ui->statusBar,SLOT(showMessage(QString,int)));
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
