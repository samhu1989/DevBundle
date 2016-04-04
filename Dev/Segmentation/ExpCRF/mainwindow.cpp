#include "mainwindow.h"
#include "ui_mainwindow.h"
#include <QMessageBox>
#include "crf2d.h"
#include <QFileDialog>
#include <QLabel>
#include "segview.h"
#include <QMdiSubWindow>
#include "densecrf.h"
MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    edit_thread_(NULL),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);
    connect(ui->actionCRF,SIGNAL(triggered(bool)),this,SLOT(start_editing()));
    connect(ui->actionCRF3D,SIGNAL(triggered(bool)),this,SLOT(start_editing()));
    connect(ui->actionLoad_Image,SIGNAL(triggered(bool)),this,SLOT(load_img()));
    connect(ui->actionLoad_Annotation,SIGNAL(triggered(bool)),this,SLOT(load_annotation()));
}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::showInMdi(QWidget* w,Qt::WindowFlags flag)
{
    w->setAttribute(Qt::WA_DeleteOnClose,true);
    QMdiSubWindow* s = ui->mdiArea->addSubWindow(w,flag);
    w->show();
    QList<QAction*> list = s->actions();
    foreach(QAction* a,list)
    {
        a->setShortcutContext(Qt::WidgetShortcut);
    }
}

void MainWindow::load_img(void)
{
    QString fileName = QFileDialog::getOpenFileName(this,
        tr("Load Input Image"),
        tr("../Dev_Data"),
        tr(";;"
        "Images Files (*.ppm,*.png,*jpg)"));
    if (fileName.isEmpty())return;
    ui->mdiArea->closeAllSubWindows();
    annotation_.clear();
    colortolabel_.clear();
    if(!input_img_.load(fileName))
    {
        QString msg = "Failed to load "+fileName+"\n";
        QMessageBox::critical(this, windowTitle(), msg);
        return;
    }
    input_img_ = input_img_.convertToFormat(QImage::Format_RGB888);
    QFileInfo info(fileName);
    input_img_.setText(tr("Path"),info.fileName());
    SegView* segview = new SegView(input_img_,annotation_,colortolabel_);
    segview->view_label(ui->actionLabel_Mask->isChecked());
    connect(ui->actionLabel_Mask,SIGNAL(triggered(bool)),segview,SLOT(view_label(bool)));
    showInMdi(segview);
    segview->update();
}

void MainWindow::load_annotation(void)
{
    QString fileName = QFileDialog::getOpenFileName(this,
        tr("Load Annotation"),
        tr("../Dev_Data"),
        tr(";;"
        "All Files (*.arma,*.ppm)"));
    if (fileName.isEmpty())return;
    QImage img(fileName);
    if(!img.isNull())
    {
        img = img.convertToFormat(QImage::Format_RGB888);
        annotation_ = DenseCRF2D::getLabelingImg(img.bits(),img.byteCount()/3,255,colortolabel_);
    }else{
        annotation_.load(fileName.toStdString());
    }
    if(annotation_.empty()){
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
        edit_thread_->setObjectName("CRF2D");
        connect(edit_thread_,SIGNAL(started()),worker,SLOT(process()));
        connect(worker,SIGNAL(end()),worker,SLOT(deleteLater()));
        connect(worker,SIGNAL(destroyed(QObject*)),edit_thread_,SLOT(quit()));
    }
    connect(edit_thread_,SIGNAL(finished()),this,SLOT(finish_editing()));
    edit_thread_->start(QThread::HighPriority);
}

void MainWindow::finish_editing()
{
    if(edit_thread_)
    {
        QString msg = edit_thread_->objectName() + "is Finished";
        QMessageBox::critical(this, windowTitle(), msg);
        edit_thread_->deleteLater();
        edit_thread_ = NULL;
    }
}
