#include "mainwindow.h"
#include "ui_mainwindow.h"
#include <QMessageBox>
#include "crf2d.h"
#include <QFileDialog>
#include <QLabel>
#include "segview.h"
#include <QMdiSubWindow>
#include "densecrf.h"
#include <OpenMesh/Tools/Utils/Timer.hh>
#include "visualizationcore.h"
#include "ncut2d.h"
MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    edit_thread_(NULL),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);
    connect(ui->actionCRF,SIGNAL(triggered(bool)),this,SLOT(start_editing()));
    connect(ui->actionCRF3D,SIGNAL(triggered(bool)),this,SLOT(start_editing()));
    connect(ui->actionNCut2D,SIGNAL(triggered(bool)),this,SLOT(start_editing()));
    connect(ui->actionLoad_Image,SIGNAL(triggered(bool)),this,SLOT(load_img()));
    connect(ui->actionLoad_Annotation,SIGNAL(triggered(bool)),this,SLOT(load_annotation()));
    connect(ui->actionLoad_Input_Mesh,SIGNAL(triggered(bool)),this,SLOT(load_mesh()));
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
    annotation_.reset(new arma::uvec());
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
    SegView* segview = new SegView(input_img_,*annotation_,colortolabel_);
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
        *annotation_ = DenseCRF2D::getLabelingImg(img.bits(),img.byteCount()/3,255,colortolabel_);
    }else{
        annotation_->load(fileName.toStdString());
    }
    if(annotation_->empty()){
        QString msg = "Failed to load "+fileName+"\n";
        QMessageBox::critical(this, windowTitle(), msg);
    }
}

void MainWindow::load_mesh(void)
{
    QString fileNames = QFileDialog::getOpenFileName(this,
        tr("Open Source file"),
        tr("../Dev_Data/"),
        tr(
        "PLY Files (*.ply);;"
        "OBJ Files (*.obj);;"
        "OFF Files (*.off);;"
        "STL Files (*.stl);;"
        "All Files (*)"));
    if (!fileNames.isEmpty())
    {
        read_mesh(fileNames);
        view_mesh();
    };
}

void MainWindow::read_mesh(const QString &filename)
{
    input_mesh_.reset(new MeshBundle<DefaultMesh>);
    ui->statusBar->showMessage(tr("Loading:")+filename,5);
    open_mesh(input_mesh_->mesh_,filename.toStdString());
    QFileInfo info(filename);
    input_mesh_->name_ = info.completeBaseName().toStdString();
}

void MainWindow::view_mesh(void)
{
    gl_timer_.stop();
    ui->mdiArea->closeAllSubWindows();
    ui->actionLabel_Mask->setChecked(false);
    MeshBundle<DefaultMesh>::Ptr& bundle_ptr = input_mesh_;
    input_mesh_view_.reset( new MeshLabelViewerWidget(this) );
    input_mesh_view_->set_mesh(bundle_ptr);
    connect(ui->actionLabel_Mask,SIGNAL(toggled(bool)),input_mesh_view_.get(),SLOT(use_custom_color(bool)));
    connect(&gl_timer_,SIGNAL(timeout()),input_mesh_view_.get(),SLOT(updateGL()));
    showInMdi((QWidget*)input_mesh_view_.get(),Qt::Widget|Qt::WindowMinMaxButtonsHint);
    gl_timer_.start(100);
}

bool MainWindow::open_mesh(DefaultMesh& mesh_,const std::string&_filename)
{
    OpenMesh::FPropHandleT< DefaultMesh::Point > fp_normal_base_;

    mesh_.request_face_normals();
    mesh_.request_face_colors();
    mesh_.request_vertex_normals();
    mesh_.request_vertex_colors();
    mesh_.request_vertex_texcoords2D();

    using namespace OpenMesh;

    IO::Options io_opt_;
    io_opt_ += OpenMesh::IO::Options::VertexColor;
    io_opt_ += OpenMesh::IO::Options::VertexNormal;
    io_opt_ += OpenMesh::IO::Options::VertexTexCoord;
    io_opt_ += OpenMesh::IO::Options::FaceColor;
    io_opt_ += OpenMesh::IO::Options::FaceNormal;
    io_opt_ += OpenMesh::IO::Options::FaceTexCoord;
    std::cout << "Loading from file '" << _filename << "'\n";
    if ( IO::read_mesh(mesh_, _filename, io_opt_ ) )
    {
      // store read option
//      io_opt_ = _opt;

      // update face and vertex normals
      if ( ! io_opt_.check( IO::Options::FaceNormal ) )
        mesh_.update_face_normals();
      else
        std::cerr << "File provides face normals\n";

      if ( ! io_opt_.check( IO::Options::VertexNormal ) )
        mesh_.update_vertex_normals();
      else
        std::cerr << "File provides vertex normals\n";


      // check for possible color information
      if ( io_opt_.check( IO::Options::VertexColor ) )
      {
        std::cerr << "File provides vertex colors\n";
      }
      else
        mesh_.release_vertex_colors();

      if ( io_opt_.check( IO::Options::FaceColor ) )
      {
        std::cerr << "File provides face colors\n";
      }
      else
        mesh_.release_face_colors();

      if ( io_opt_.check( IO::Options::VertexTexCoord ) )
        std::cerr << "File provides texture coordinates\n";

      // info
      std::clog << mesh_.n_vertices() << " vertices, "
            << mesh_.n_edges()    << " edge, "
            << mesh_.n_faces()    << " faces\n";

      // base point for displaying face normals
      {
        OpenMesh::Utils::Timer t;
        t.start();
        mesh_.add_property( fp_normal_base_ );
        DefaultMesh::FaceIter f_it = mesh_.faces_begin();
        DefaultMesh::FaceVertexIter fv_it;
        for (;f_it != mesh_.faces_end(); ++f_it)
        {
          DefaultMesh::Point v(0,0,0);
          for( fv_it=mesh_.fv_iter(*f_it); fv_it.is_valid(); ++fv_it)
            v += OpenMesh::vector_cast<DefaultMesh::Normal>(mesh_.point(*fv_it));
          v *= 1.0f/3.0f;
          mesh_.property( fp_normal_base_, *f_it ) = v;
        }
        t.stop();
        std::clog << "Computed base point for displaying face normals ["
                  << t.as_string() << "]" << std::endl;
      }
      return true;
    }
    return false;
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
    if(annotation_->empty())
    {
        QString msg = "Please Load Annotation First\n";
        QMessageBox::critical(this, windowTitle(), msg);
        return;
    }
    QAction* edit = qobject_cast<QAction*>(sender());
    if(edit==ui->actionCRF)
    {
        edit_thread_ = new QThread();
        CRF2D* worker = new CRF2D(input_img_,*annotation_);
        connect(worker,SIGNAL(message(QString,int)),ui->statusBar,SLOT(showMessage(QString,int)));
        worker->moveToThread(edit_thread_);
        edit_thread_->setObjectName("CRF2D");
        connect(edit_thread_,SIGNAL(started()),worker,SLOT(process()));
        connect(worker,SIGNAL(end()),worker,SLOT(deleteLater()));
        connect(worker,SIGNAL(destroyed(QObject*)),edit_thread_,SLOT(quit()));
    }
    if(edit==ui->actionNCut2D)
    {
        edit_thread_ = new QThread();
        NCut2D* worker = new NCut2D(input_img_,*annotation_);
        connect(worker,SIGNAL(message(QString,int)),ui->statusBar,SLOT(showMessage(QString,int)));
        worker->moveToThread(edit_thread_);
        edit_thread_->setObjectName("NCut2D");
        connect(edit_thread_,SIGNAL(started()),worker,SLOT(process()));
        connect(worker,SIGNAL(end()),edit_thread_,SLOT(quit()));
        connect(worker,SIGNAL(end()),worker,SLOT(deleteLater()));
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
