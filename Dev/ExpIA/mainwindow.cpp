#include "mainwindow.h"
#include "ui_mainwindow.h"
#include <QFileDialog>
#include <QFileInfo>
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Tools/Utils/Timer.hh>
#include "MeshPairViewerWidget.h"
#include "regiongrowthread.h"
#include <vector>
MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow),
    edit_thread_(NULL),config_(new Config("./Default.config"))
{
    ui->setupUi(this);

    if(!config_->has("Configure"))config_->reload("../Default.config");

    connect(ui->actionConfigure,SIGNAL(triggered(bool)),this,SLOT(configure()));
    connect(ui->actionOpen_Inputs,SIGNAL(triggered(bool)),this,SLOT(open_inputs()));
    connect(ui->actionRegionGrow,SIGNAL(triggered(bool)),this,SLOT(start_editing()));
    io_opt_ += OpenMesh::IO::Options::VertexColor;
    io_opt_ += OpenMesh::IO::Options::VertexNormal;
    io_opt_ += OpenMesh::IO::Options::VertexTexCoord;
    io_opt_ += OpenMesh::IO::Options::FaceColor;
    io_opt_ += OpenMesh::IO::Options::FaceNormal;
    io_opt_ += OpenMesh::IO::Options::FaceTexCoord;
}

void MainWindow::configure()
{
    QString fileName = QFileDialog::getOpenFileName(this,
        tr("Confugre"),
        tr("./"),
        tr("Configure (*.config);;"
        "All Files (*)"));
    if (!fileName.isEmpty())
    {
        config_->reload(fileName.toStdString());
    }
}

void MainWindow::open_inputs()
{
    QStringList fileNames = QFileDialog::getOpenFileNames(this,
        tr("Open Source file"),
        tr("../Dev_Data/"),
        tr("OBJ Files (*.obj);;"
        "OFF Files (*.off);;"
        "STL Files (*.stl);;"
        "PLY Files (*.ply);;"
        "All Files (*)"));
    if (!fileNames.isEmpty())
    {
        open_inputs(fileNames);
        view_inputs();
    }
}

void MainWindow::open_inputs(QStringList&fileNames)
{
    inputs_.clear();
    foreach(QString fname,fileNames)
    {
        ui->statusBar->showMessage(tr("Loading:")+fname,5);
        inputs_.push_back(std::make_shared<MeshBundle<DefaultMesh>>());
        open_mesh(inputs_.back()->mesh_,fname.toStdString());
        QFileInfo info(fname);
        inputs_.back()->name_ = info.completeBaseName().toStdString();
        QApplication::processEvents();
    }
}

bool MainWindow::open_mesh(DefaultMesh& mesh_,const std::string&_filename)
{
    OpenMesh::FPropHandleT< DefaultMesh::Point > fp_normal_base_;

    mesh_.request_face_normals();
    mesh_.request_face_colors();
    mesh_.request_vertex_normals();
    mesh_.request_vertex_colors();
    mesh_.request_vertex_texcoords2D();

    IO::Options _opt = io_opt_;
    std::cout << "Loading from file '" << _filename << "'\n";
    if ( IO::read_mesh(mesh_, _filename, _opt ) )
    {
      // store read option
      io_opt_ = _opt;

      // update face and vertex normals
      if ( ! io_opt_.check( IO::Options::FaceNormal ) )
        mesh_.update_face_normals();
      else
        std::cout << "File provides face normals\n";

      if ( ! io_opt_.check( IO::Options::VertexNormal ) )
        mesh_.update_vertex_normals();
      else
        std::cout << "File provides vertex normals\n";


      // check for possible color information
      if ( io_opt_.check( IO::Options::VertexColor ) )
      {
        std::cout << "File provides vertex colors\n";
      }
      else
        mesh_.release_vertex_colors();

      if ( io_opt_.check( IO::Options::FaceColor ) )
      {
        std::cout << "File provides face colors\n";
      }
      else
        mesh_.release_face_colors();

      if ( io_opt_.check( IO::Options::VertexTexCoord ) )
        std::cout << "File provides texture coordinates\n";

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

void MainWindow::view_inputs()
{
    ui->mdiArea->closeAllSubWindows();
    mesh_views_.clear();
    ui->actionCustom_Color->setChecked(false);
    std::vector<MeshBundle<DefaultMesh>::Ptr>::iterator iter;
    for(iter=inputs_.begin();iter!=inputs_.end();++iter)
    {
        MeshBundle<DefaultMesh>::Ptr& bundle_ptr = *iter;
        MeshPairViewerWidget* widget = new MeshPairViewerWidget(this);
        widget->setMinimumSize(300,200);
        widget->first_ptr() = bundle_ptr;
        widget->set_center_at_mesh(bundle_ptr->mesh_);
        widget->setWindowTitle(QString::fromStdString(bundle_ptr->name_));
        connect(ui->actionCustom_Color,SIGNAL(toggled(bool)),widget,SLOT(use_custom_color(bool)));
        mesh_views_.push_back(WidgetPtr(widget));
        showInMdi((QWidget*)widget);
    }
}

void MainWindow::showInMdi(QWidget* w)
{
    w->setAttribute(Qt::WA_DeleteOnClose,true);
    ui->mdiArea->addSubWindow(w);
    connect(w,SIGNAL(destroyed()),this,SLOT(removeView()));
    w->show();
}

void MainWindow::removeView()
{
    WidgetPtr w = qobject_cast<WidgetPtr>(sender());
    if(w)
    {
        if(mesh_views_.back()==w)mesh_views_.pop_back();
        else{
            std::vector<WidgetPtr>::iterator viter;
            for(viter=mesh_views_.begin();viter!=mesh_views_.end();++viter)
            {
                if(*viter==w)
                {
                    *viter = mesh_views_.back();
                    mesh_views_.pop_back();
                    break;
                }
            }
        }
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
    QAction* edit = qobject_cast<QAction*>(sender());
    if(edit==ui->actionRegionGrow)
    {
        RegionGrowThread* th = new RegionGrowThread(inputs_,labels_);
        if(!th->configure(config_)){
            QString msg = "Missing Some Configure\n";
            QMessageBox::critical(this, windowTitle(), msg);
            return;
        }
        connect(th,SIGNAL(message(QString,int)),ui->statusBar,SLOT(showMessage(QString,int)));
        edit_thread_ = th;
    }
    connect(edit_thread_,SIGNAL(finished()),this,SLOT(finish_editing()));
    if(edit_thread_)edit_thread_->start(QThread::HighestPriority);
}

void MainWindow::finish_editing()
{
    QString msg = edit_thread_->objectName() + " is Finished";
    if(edit_thread_)
    {
        while(edit_thread_->isRunning())
        {
            edit_thread_->terminate();
            QApplication::processEvents();
        }
        edit_thread_->deleteLater();
        edit_thread_ = NULL;
    }
    QMessageBox::information( this, windowTitle(), msg);
}

MainWindow::~MainWindow()
{
    delete ui;
}
