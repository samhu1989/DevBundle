#include "mainwindow.h"
#include "ui_mainwindow.h"
#include <QFileDialog>
#include <QFileInfo>
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Tools/Utils/Timer.hh>
#include "MeshPairViewerWidget.h"
#include "regiongrowthread.h"
#include "unifylabelthread.h"
#include "unifylabelmannual.h"
#include "updateobjectmodel.h"
#include "labspace.h"
#include <vector>
#include <QMdiSubWindow>
#include <QDebug>
#include "supervoxelthread.h"
#include "graphcutthread.h"
#include "objectview.h"
#include "globalalign.h"
#include "extractbackground.h"
#include "inpatchgraphcut.h"
#include <typeinfo>
#include <fstream>
#include <strstream>
MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow),
    edit_thread_(NULL),config_(new Config("./Default.config"))
{
    ui->setupUi(this);

    gl_timer.stop();
    if(!config_->has("Configure"))config_->reload("../Default.config");

    connect(ui->actionConfigure,SIGNAL(triggered(bool)),this,SLOT(configure()));
    connect(ui->actionOpen_Inputs,SIGNAL(triggered(bool)),this,SLOT(open_inputs()));
    connect(ui->actionSave_Aligned,SIGNAL(triggered(bool)),this,SLOT(save_aligned()));
    connect(ui->actionSave_Supervoxels,SIGNAL(triggered(bool)),this,SLOT(save_supervoxels()));
    connect(ui->actionLoad_Supervoxels,SIGNAL(triggered(bool)),this,SLOT(load_supervoxels()));
    connect(ui->actionSave_Segments,SIGNAL(triggered(bool)),this,SLOT(save_labels()));
    connect(ui->actionLoad_Segments,SIGNAL(triggered(bool)),this,SLOT(load_labels()));
    connect(ui->actionSave_Object_Model,SIGNAL(triggered(bool)),this,SLOT(save_objects()));
    connect(ui->actionLoad_Objects,SIGNAL(triggered(bool)),this,SLOT(load_objects()));
    connect(ui->actionSave_Scenes,SIGNAL(triggered(bool)),this,SLOT(save_scenes()));

    connect(ui->actionGlobal_Align,SIGNAL(triggered(bool)),this,SLOT(start_editing()));
    connect(ui->actionExtract_Background,SIGNAL(triggered(bool)),this,SLOT(start_editing()));
    connect(ui->actionRemove_Zero_Label,SIGNAL(triggered(bool)),this,SLOT(remove_zero_label()));
    connect(ui->actionSupervoxel,SIGNAL(triggered(bool)),this,SLOT(start_editing()));
    connect(ui->actionRegionGrow,SIGNAL(triggered(bool)),this,SLOT(start_editing()));
    connect(ui->actionAutomatically,SIGNAL(triggered(bool)),this,SLOT(start_editing()));
    connect(ui->actionMannually,SIGNAL(triggered(bool)),this,SLOT(start_editing()));
    connect(ui->actionUpdateObjModel,SIGNAL(triggered(bool)),this,SLOT(start_editing()));
    connect(ui->actionGlobal_Graph_Cut,SIGNAL(triggered(bool)),this,SLOT(start_editing()));
    connect(ui->actionIn_Patch_Graph_Cut,SIGNAL(triggered(bool)),this,SLOT(start_editing()));

    connect(ui->actionLab_Color_Space,SIGNAL(triggered(bool)),this,SLOT(showLab()));
    connect(ui->actionObject_View,SIGNAL(triggered(bool)),this,SLOT(viewObj()));
    connect(ui->actionSupervoxel_Color,SIGNAL(triggered(bool)),this,SLOT(showSVColor()));

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
        tr(
        "PLY Files (*.ply);;"
        "OBJ Files (*.obj);;"
        "OFF Files (*.off);;"
        "STL Files (*.stl);;"
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
    labels_.clear();
    foreach(QString fname,fileNames)
    {
        ui->statusBar->showMessage(tr("Loading:")+fname,5);
        inputs_.push_back(std::make_shared<MeshBundle<DefaultMesh>>());
        open_mesh(inputs_.back()->mesh_,fname.toStdString());
        QFileInfo info(fname);
        inputs_.back()->name_ = info.completeBaseName().toStdString();
        labels_.emplace_back(inputs_.back()->mesh_.n_vertices(),arma::fill::zeros);
        QApplication::processEvents();
    }
}

void MainWindow::save_aligned()
{
    QString dirName = QFileDialog::getExistingDirectory(
                this,
                tr("Export Aligned Inputs"),
                tr("../Dev_Data/")
                );
    if(dirName.isEmpty())return;
    OpenMesh::IO::Options opt;
    opt+=OpenMesh::IO::Options::Binary;
    opt+=OpenMesh::IO::Options::VertexColor;
    opt+=OpenMesh::IO::Options::VertexNormal;
    std::string path = dirName.toStdString();
    MeshBundle<DefaultMesh>::PtrList::iterator iter;
    for(iter=inputs_.begin();iter!=inputs_.end();++iter)
    {
        MeshBundle<DefaultMesh>::Ptr ptr = *iter;
        if(!OpenMesh::IO::write_mesh(ptr->mesh_,path+"/"+ptr->name_+"_aligned.ply",opt,13)){
            std::cerr<<"can't save to:"<<path+"/"+ptr->name_+"_aligned.ply"<<std::endl;
            return;
        }
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
    gl_timer.stop();
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
        connect(&gl_timer,SIGNAL(timeout()),widget,SLOT(updateGL()));
        mesh_views_.push_back(WidgetPtr(widget));
        showInMdi((QWidget*)widget,Qt::Widget|Qt::WindowMinMaxButtonsHint);
    }
    gl_timer.start(100);
}

void MainWindow::save_labels()
{
    QString dirName = QFileDialog::getExistingDirectory(
                this,
                tr("Save Labels"),
                tr("../Dev_Data/")
                );
    if(dirName.isEmpty())return;
    std::vector<arma::uvec>::iterator iter;
    std::vector<MeshBundle<DefaultMesh>::Ptr>::iterator miter;
    QDir dir;
    dir.setPath(dirName);
    miter = inputs_.begin();
    for(iter=labels_.begin();iter!=labels_.end();++iter)
    {
        QString filepath = dir.absoluteFilePath(
                    QString::fromStdString((*miter)->name_+".label.arma")
                    );
        if(!iter->save(filepath.toStdString(),arma::arma_binary))
        {
            QString msg = "Failed to Save "+filepath+"\n";
            QMessageBox::critical(this, windowTitle(), msg);
            return;
        }
        if(miter==inputs_.end())break;
        ++miter;
    }
}

void MainWindow::load_labels()
{
    if(inputs_.empty())
    {
        QString msg = "Please Load Inputs First\n";
        QMessageBox::critical(this, windowTitle(), msg);
    }
    QString dirName = QFileDialog::getExistingDirectory(
                this,
                tr("Load Labels"),
                tr("../Dev_Data/")
                );
    if(dirName.isEmpty())return;
    std::vector<arma::uvec>::iterator iter;
    std::vector<MeshBundle<DefaultMesh>::Ptr>::iterator miter;
    QDir dir;
    dir.setPath(dirName);
    miter = inputs_.begin();
    if(labels_.size()!=inputs_.size())
    {
        labels_.resize(inputs_.size());
    }
    for(iter=labels_.begin();iter!=labels_.end();++iter)
    {
        arma::uvec label;
        QString filepath = dir.absoluteFilePath(
                    QString::fromStdString((*miter)->name_+".label.arma")
                    );
        if(!label.load(filepath.toStdString()))
        {
            QString msg = "Failed to Load "+filepath+"\n";
            QMessageBox::critical(this, windowTitle(), msg);
            return;
        }
        if( label.size() != (*miter)->mesh_.n_vertices() )
        {
            QString msg = "Some size of this set of labels doesn't match the inputs\n";
            QMessageBox::critical(this, windowTitle(), msg);
            return;
        }
        *iter = label;
        (*miter)->custom_color_.fromlabel(*iter);
        if(miter==inputs_.end())break;
        ++miter;
    }
}

void MainWindow::save_objects()
{
    QString dirName = QFileDialog::getExistingDirectory(
                this,
                tr("Save Objects"),
                tr("../Dev_Data/")
                );
    if(dirName.isEmpty())return;
    std::vector<ObjModel::Ptr>::iterator iter;
    QDir dir;
    dir.setPath(dirName);
    size_t index = 0;
    for(iter=objects_.begin();iter!=objects_.end();++iter)
    {
        ObjModel::Ptr &ptr = *iter;
        if( !ptr || 0 == ptr.use_count()){
            std::cerr<<"Not all Objects Generated"<<std::endl;
            break;
        }else if( 0 == ptr->GeoM_->mesh_.n_vertices() )
        {
            std::cerr<<"Empty Objects In List"<<std::endl;
            break;
        }
        QString path;
        path = path.sprintf("GeoObj%u",index);
        QDir odir;
        odir.setPath(dir.absoluteFilePath(path));
        if(!odir.exists())
        {
            if(!dir.mkdir(path))
            {
                QString msg = "Failed to make "+dir.absoluteFilePath(path)+"\n";
                QMessageBox::critical(this, windowTitle(), msg);
                return;
            }
        }
        QString filepath = dir.absoluteFilePath(path);
        if(!ptr->save(filepath.toStdString()))
        {
            QString msg = "Failed to Save to"+filepath+"\n";
            QMessageBox::critical(this, windowTitle(), msg);
            return;
        }
        ++index;
    }
}

void MainWindow::load_objects()
{
    QString dirName = QFileDialog::getExistingDirectory(
                this,
                tr("Load Objects"),
                tr("../Dev_Data/")
                );
    if(dirName.isEmpty())return;
    std::vector<ObjModel::Ptr>::iterator iter;
    QDir dir;
    dir.setPath(dirName);
    QString path;
    path = "GeoObj0";
    QDir odir;
    odir.setPath(dir.absoluteFilePath(path));
    objects_.clear();
    size_t index = 0;
    while(odir.exists())
    {
        objects_.emplace_back(new ObjModel());
        ObjModel::Ptr &ptr = objects_.back();
        if(!ptr->load(odir.absolutePath().toStdString()))
        {
            QString msg = "Failed to Load From"+odir.absolutePath()+"\n";
            QMessageBox::critical(this, windowTitle(), msg);
            return;
        }
        ++index;
        path = path.sprintf("GeoObj%u",index);
        QString filepath = dir.absoluteFilePath(path);
        odir.setPath(filepath);
    }
    if(!labels_.empty())
    {
        size_t max_label_ = 0;
        std::vector<arma::uvec>::iterator iter;
        for(iter=labels_.begin();iter!=labels_.end();++iter)
        {
            arma::uword max = arma::max(*iter);
            if(max>max_label_)max_label_=max;
        }
        if(max_label_!=objects_.size())
        {
            QString msg;
            msg = msg.sprintf("%u objects loaded \n while labels show %u objects",objects_.size(),max_label_);
            QMessageBox::warning(this, windowTitle(), msg);
            return;
        }
    }
}

void MainWindow::save_supervoxels()
{
    QString dirName = QFileDialog::getExistingDirectory(
                this,
                tr("Save Supervoxels"),
                tr("../Dev_Data/")
                );
    if(dirName.isEmpty())return;
    std::vector<MeshBundle<DefaultMesh>::Ptr>::iterator iter;
    QDir dir;
    dir.setPath(dirName);
    for(iter=inputs_.begin();iter!=inputs_.end();++iter)
    {
        MeshBundle<DefaultMesh>::Ptr &ptr = *iter;
        if( !ptr || 0 == ptr.use_count()){
            std::cerr<<"Empty Input ?"<<std::endl;
            break;
        }else if( 0 == ptr->graph_.voxel_neighbors.size() ||
                  0 == ptr->graph_.voxel_label.size()   ||
                  0 == ptr->graph_.voxel_centers.size()
                  )
        {
            std::cerr<<"Empty Supervoxel In List"<<std::endl;
            break;
        }
        QString path;
        path = "supervoxel_"+QString::fromStdString(ptr->name_);
        QDir odir;
        odir.setPath(dir.absoluteFilePath(path));
        if(!odir.exists())
        {
            if(!dir.mkdir(path))
            {
                QString msg = "Failed to make "+dir.absoluteFilePath(path)+"\n";
                QMessageBox::critical(this, windowTitle(), msg);
                return;
            }
        }
        QString filepath = dir.absoluteFilePath(path);
        if(!ptr->graph_.save(filepath.toStdString()))
        {
            QString msg = "Failed to Save to"+filepath+"\n";
            QMessageBox::critical(this, windowTitle(), msg);
            return;
        }
    }
}

void MainWindow::load_supervoxels()
{
    QString dirName = QFileDialog::getExistingDirectory(
                this,
                tr("Load Supervoxels"),
                tr("../Dev_Data/")
                );
    if(dirName.isEmpty())return;
    std::vector<MeshBundle<DefaultMesh>::Ptr>::iterator iter;
    QDir dir;
    dir.setPath(dirName);
    for(iter=inputs_.begin();iter!=inputs_.end();++iter)
    {
        MeshBundle<DefaultMesh>::Ptr &ptr = *iter;
        if( !ptr || 0 == ptr.use_count()){
            std::cerr<<"Empty Input ?"<<std::endl;
            break;
        }
        QString path;
        path = "supervoxel_"+QString::fromStdString(ptr->name_);
        QDir odir;
        odir.setPath(dir.absoluteFilePath(path));
        if(!odir.exists())
        {
            QString msg = "Failed to find "+dir.absoluteFilePath(path)+"\n";
            QMessageBox::critical(this, windowTitle(), msg);
            return;
        }
        QString filepath = dir.absoluteFilePath(path);
        if(!ptr->graph_.load(filepath.toStdString()))
        {
            QString msg = "Failed to Load from"+filepath+"\n";
            QMessageBox::critical(this, windowTitle(), msg);
            return;
        }
        ptr->custom_color_.fromlabel(ptr->graph_.voxel_label);
    }
}

void MainWindow::save_scenes()
{
    QString dirName = QFileDialog::getExistingDirectory(
                this,
                tr("Save Scenes"),
                tr("../Dev_Data/")
                );
    if(dirName.isEmpty())return;
    if(inputs_.empty())return;
    if(objects_.empty())return;
    if(!config_->has("O_model_suffix")){
        QMessageBox::information( this, windowTitle(), tr("Missing Config \"O_model_suffix\" "));
        return;
    }
    QDir dir;
    dir.setPath(dirName+"/"+QString::fromStdString("layouts/"));
    if(!dir.exists())
    {
        if(!dir.mkdir(dirName+"/"+QString::fromStdString("layouts/")))
        {
            QString msg = "Failed to make "+dir.absolutePath()+"\n";
            QMessageBox::critical(this, windowTitle(), msg);
            return;
        }
    }
    //object layouts
    save_object_layout(dirName.toStdString());
    //scene layouts
    save_scene_layout(dirName.toStdString());
    //scene model
    dir.setPath(dirName+"/"+QString::fromStdString("models/"));
    if(!dir.exists())
    {
        if(!dir.mkdir(dirName+"/"+QString::fromStdString("models/")))
        {
            QString msg = "Failed to make "+dir.absolutePath()+"\n";
            QMessageBox::critical(this, windowTitle(), msg);
            return;
        }
    }
    save_scene_model(dirName.toStdString());
}

void MainWindow::save_object_layout(const std::string&path)
{
    std::fstream objfile;
    std::stringstream stream("");
    std::stringstream layout_stream("");
    std::string objfilename;
    size_t index = 1;
    std::vector<ObjModel::Ptr>::iterator iter;
    for(iter=objects_.begin();iter!=objects_.end();++iter)
    {
        ObjModel::Ptr obj_ptr = *iter;
        stream.clear();
        stream<<path<<"/layouts/obj_box"<<index<<".obj";
        stream>>objfilename;
        objfile.open(objfilename,objfile.out);
        if(!objfile.is_open())
        {
            std::cerr<<"Failed to Open:"<<objfilename<<std::endl;
        }
        objfile<<"#obj "<<index<<std::endl;
        arma::fmat layout_mat;
        if(obj_ptr->fullLayout(layout_mat,-1));
        {
            layout_stream.clear();
            size_t cnt = 1;
            for(size_t c=0;c<layout_mat.n_cols;++c)
            {
                layout_stream<<"v "<<layout_mat(0,c)
                            <<" "<<layout_mat(1,c)
                           <<" "<<layout_mat(2,c)
                          <<std::endl;
            }
            layout_stream<<"f "<<cnt+0<<" "<<cnt+1<<" "<<cnt+2<<" "<<cnt+3<<std::endl;
            layout_stream<<"f "<<cnt+4<<" "<<cnt+5<<" "<<cnt+6<<" "<<cnt+7<<std::endl;
            layout_stream<<"f "<<cnt+0<<" "<<cnt+4<<" "<<cnt+7<<" "<<cnt+3<<std::endl;
            layout_stream<<"f "<<cnt+0<<" "<<cnt+4<<" "<<cnt+5<<" "<<cnt+1<<std::endl;
            layout_stream<<"f "<<cnt+1<<" "<<cnt+5<<" "<<cnt+6<<" "<<cnt+2<<std::endl;
            layout_stream<<"f "<<cnt+2<<" "<<cnt+6<<" "<<cnt+7<<" "<<cnt+3<<std::endl;
            objfile << layout_stream.str();
            objfile.close();
        }
        ++index;
    }
}

void MainWindow::save_scene_layout(const std::string&path)
{
    std::fstream objfile;
    std::stringstream stream("");
    std::string objfilename;
    for(size_t fIndex = 0 ; fIndex < inputs_.size() ; ++fIndex )
    {
        stream.clear();
        stream.str("");
        stream<<path<<"/layouts/scene_"<<fIndex<<".obj";
        stream>>objfilename;
        objfile.open(objfilename,objfile.out);
        if(!objfile.is_open())
        {
            std::cerr<<"Failed to Open:"<<objfilename<<std::endl;
        }
        size_t cnt = 1;
        std::vector<ObjModel::Ptr>::iterator iter;
        size_t oIndex = 1;
        for(iter=objects_.begin();iter!=objects_.end();++iter)
        {
            objfile<<"#obj "<<oIndex<<std::endl;
            arma::fmat layout_mat;
            if((*iter)->fullLayout(layout_mat,fIndex));
            {
                for(size_t c=0;c<layout_mat.n_cols;++c)
                {
                    objfile<<"v "<<layout_mat(0,c)
                          <<" "<<layout_mat(1,c)
                         <<" "<<layout_mat(2,c)
                        <<std::endl;
                }
                objfile<<"f "<<cnt+0<<" "<<cnt+1<<" "<<cnt+2<<" "<<cnt+3<<std::endl;
                objfile<<"f "<<cnt+4<<" "<<cnt+5<<" "<<cnt+6<<" "<<cnt+7<<std::endl;
                objfile<<"f "<<cnt+0<<" "<<cnt+4<<" "<<cnt+7<<" "<<cnt+3<<std::endl;
                objfile<<"f "<<cnt+0<<" "<<cnt+4<<" "<<cnt+5<<" "<<cnt+1<<std::endl;
                objfile<<"f "<<cnt+1<<" "<<cnt+5<<" "<<cnt+6<<" "<<cnt+2<<std::endl;
                objfile<<"f "<<cnt+2<<" "<<cnt+6<<" "<<cnt+7<<" "<<cnt+3<<std::endl;
                cnt+=8;
            }
            ++oIndex;
        }
        objfile.close();
    }
}

void MainWindow::save_scene_model(const std::string&path)
{
    size_t index = 1;
    std::stringstream obj_stream("");
    OpenMesh::IO::Options opt;
    opt+=OpenMesh::IO::Options::Binary;
    opt+=OpenMesh::IO::Options::VertexColor;
    opt+=OpenMesh::IO::Options::VertexNormal;
    std::vector<ObjModel::Ptr>::iterator iter;
    for(iter=objects_.begin();iter!=objects_.end();++iter)
    {
        ObjModel::Ptr obj_ptr = *iter;
        DefaultMesh m;
        obj_ptr->fullModel(m,-1);
        obj_stream.clear();
        obj_stream.str("");
        obj_stream<<path+"/models/obj"<<index<<config_->getString("O_model_suffix");
        if(!OpenMesh::IO::write_mesh(m,obj_stream.str(),opt,10)){
            std::cerr<<"can't save to:"<<obj_stream.str()<<std::endl;
            return;
        }
        ++index;
    }
}

void MainWindow::showInMdi(QWidget* w,Qt::WindowFlags flag)
{
    w->setAttribute(Qt::WA_DeleteOnClose,true);
    QMdiSubWindow* s = ui->mdiArea->addSubWindow(w,flag);
    connect(w,SIGNAL(destroyed()),this,SLOT(removeView()));
    w->show();
    QList<QAction*> list = s->actions();
    foreach(QAction* a,list)
    {
        a->setShortcutContext(Qt::WidgetShortcut);
    }
}

void MainWindow::closeInMdi(QWidget *w)
{
    ui->mdiArea->setActiveSubWindow((QMdiSubWindow*)w);
    ui->mdiArea->closeActiveSubWindow();
}

void MainWindow::showBox(int index,MeshBundle<DefaultMesh>::Ptr ptr)
{
    if(index>=mesh_views_.size())return;
    MeshPairViewerWidget* w = qobject_cast<MeshPairViewerWidget*>(mesh_views_[index]);
    w->second_ptr() = ptr;
}

void MainWindow::removeView()
{
    WidgetPtr w = qobject_cast<WidgetPtr>(sender());
    if(w)
    {
        if(mesh_views_.empty())return;
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
    if( inputs_.empty() || labels_.empty() )
    {
        QString msg = "No Inputs for Editing\n";
        QMessageBox::critical(this, windowTitle(), msg);
        return;
    }
    QAction* edit = qobject_cast<QAction*>(sender());
    if(edit==ui->actionSupervoxel)
    {
        SupervoxelThread* th = new SupervoxelThread(inputs_);
        if(!th->configure(config_)){
            th->deleteLater();
            QString msg = "Missing Some Configure\n";
            QMessageBox::critical(this, windowTitle(), msg);
            return;
        }
        connect(th,SIGNAL(message(QString,int)),ui->statusBar,SLOT(showMessage(QString,int)));
        edit_thread_ = th;
    }
    if(edit==ui->actionRegionGrow)
    {
        RegionGrowThread* th = new RegionGrowThread(inputs_,labels_);
        if(!th->configure(config_)){
            th->deleteLater();
            QString msg = "Missing Some Configure\n";
            QMessageBox::critical(this, windowTitle(), msg);
            return;
        }
        connect(th,SIGNAL(message(QString,int)),ui->statusBar,SLOT(showMessage(QString,int)));
        edit_thread_ = th;
    }
    if(edit==ui->actionAutomatically)
    {
        UnifyLabelThread* th = new UnifyLabelThread(inputs_,labels_,feature_base_);
        if(!th->configure(config_)){
            th->deleteLater();
            QString msg = "Missing Some Configure\n";
            QMessageBox::critical(this, windowTitle(), msg);
            return;
        }
        connect(th,SIGNAL(message(QString,int)),ui->statusBar,SLOT(showMessage(QString,int)));
        edit_thread_ = th;
    }
    if(edit==ui->actionMannually)
    {
        UnifyLabelMannual* w = new UnifyLabelMannual(inputs_,labels_);
        if(!w->configure(config_)){
            QString msg = "You probably should do regiongrow first\n";
            QMessageBox::critical(this, windowTitle(), msg);
            w->deleteLater();
            return;
        }
        connect(w,SIGNAL(message(QString,int)),ui->statusBar,SLOT(showMessage(QString,int)));
        showInMdi((QWidget*)w);
        w->initLater();
    }
    if(edit==ui->actionUpdateObjModel)
    {
        UpdateObjectModel* w = new UpdateObjectModel(
                    inputs_,
                    labels_,
                    objects_
                    );
        if(!w->configure(config_)){
            QString msg = "You probably should do unify label first\n";
            QMessageBox::critical(this, windowTitle(), msg);
            w->deleteLater();
            return;
        }
        connect(w,SIGNAL(message(QString,int)),ui->statusBar,SLOT(showMessage(QString,int)));
        connect(w,SIGNAL(show_layout(int,MeshBundle<DefaultMesh>::Ptr)),this,SLOT(showBox(int,MeshBundle<DefaultMesh>::Ptr)));
        w->setAttribute(Qt::WA_DeleteOnClose,true);
        QMdiSubWindow* s = ui->mdiArea->addSubWindow(w,Qt::Widget|Qt::WindowMinMaxButtonsHint);
        connect(w,SIGNAL(closeInMdi(QWidget*)),this,SLOT(closeInMdi(QWidget*)));
        s->show();
        w->startLater();
    }
    if(edit==ui->actionGlobal_Align)
    {
        GlobalAlign* w = new GlobalAlign(inputs_);
        if(!w->configure(config_)){
            QString msg = "You probably should open inputs first\n";
            QMessageBox::critical(this, windowTitle(), msg);
            w->deleteLater();
            return;
        }
        connect(w,SIGNAL(message(QString,int)),ui->statusBar,SLOT(showMessage(QString,int)));
        w->setAttribute(Qt::WA_DeleteOnClose,true);
        QMdiSubWindow* s = ui->mdiArea->addSubWindow(w);
        s->show();
    }
    if(edit==ui->actionExtract_Background)
    {
        ExtractBackground* w = new ExtractBackground(mesh_views_,gl_timer,labels_);
        if(!w->configure(config_)){
            QString msg = "You probably should open inputs first\n";
            QMessageBox::critical(this, windowTitle(), msg);
            w->deleteLater();
            return;
        }
        connect(w,SIGNAL(message(QString,int)),ui->statusBar,SLOT(showMessage(QString,int)));
        w->setAttribute(Qt::WA_DeleteOnClose,true);
        QMdiSubWindow* s = ui->mdiArea->addSubWindow(w);
        s->show();
    }
    if(edit==ui->actionGlobal_Graph_Cut)
    {
        GraphCutThread* th = new GraphCutThread(inputs_,objects_,labels_);
        if(!th->configure(config_)){
            th->deleteLater();
            QString msg = "Missing Some Inputs or configure\n";
            QMessageBox::critical(this, windowTitle(), msg);
            return;
        }
        connect(th,SIGNAL(message(QString,int)),ui->statusBar,SLOT(showMessage(QString,int)));
        connect(th,SIGNAL(sendMatch(int,MeshBundle<DefaultMesh>::Ptr)),this,SLOT(showBox(int,MeshBundle<DefaultMesh>::Ptr)),Qt::DirectConnection);
        edit_thread_ = th;
    }
    if(edit==ui->actionIn_Patch_Graph_Cut)
    {
        InPatchGraphCut* th = new InPatchGraphCut(inputs_,objects_,labels_);
        if(!th->configure(config_)){
            th->deleteLater();
            QString msg = "Missing Some Inputs or configure\n";
            QMessageBox::critical(this, windowTitle(), msg);
            return;
        }
        connect(th,SIGNAL(message(QString,int)),ui->statusBar,SLOT(showMessage(QString,int)));
        connect(th,SIGNAL(sendMatch(int,MeshBundle<DefaultMesh>::Ptr)),this,SLOT(showBox(int,MeshBundle<DefaultMesh>::Ptr)),Qt::DirectConnection);
        edit_thread_ = th;
    }
    if(edit_thread_){
        connect(edit_thread_,SIGNAL(finished()),this,SLOT(finish_editing()));
        edit_thread_->start(QThread::HighestPriority);
    }
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

void MainWindow::remove_zero_label()
{
    MeshBundle<DefaultMesh>::PtrList::iterator iter;
    std::vector<arma::uvec>::iterator liter = labels_.begin();
    size_t mindex = 0;
    for(iter=inputs_.begin();iter!=inputs_.end();++iter)
    {
        MeshPairViewerWidget* w =  (MeshPairViewerWidget*)mesh_views_[mindex];
        disconnect(&gl_timer,SIGNAL(timeout()),w,SLOT(updateGL()));
        DefaultMesh mesh;
        mesh.request_vertex_colors();
        mesh.request_vertex_normals();
        MeshBundle<DefaultMesh>::Ptr ptr = (*iter);
        arma::uvec& label = (*liter);
        long long index = -1;
        for (DefaultMesh::VertexIter v_it = ptr->mesh_.vertices_begin();v_it != ptr->mesh_.vertices_end();++v_it)
        {
            ++index;
            if(0==label(index))continue;
            DefaultMesh::VertexHandle new_it = mesh.add_vertex(ptr->mesh_.point(*v_it));
            mesh.set_color(new_it,ptr->mesh_.color(*v_it));
            mesh.set_normal(new_it,ptr->mesh_.normal(*v_it));
        }
        ptr->mesh_ = mesh;
        arma::uvec label_indices = arma::find( label != 0 );
        label = label(label_indices);
        ptr->custom_color_.fromlabel(label);
        connect(&gl_timer,SIGNAL(timeout()),w,SLOT(updateGL()));
        ++liter;
        ++mindex;
        if(liter==labels_.end())break;
        QApplication::processEvents();
    }
}

void MainWindow::showLab()
{
    LabSpace* v = new LabSpace();
    v->setAttribute(Qt::WA_DeleteOnClose,true);
    QMdiSubWindow* s = ui->mdiArea->addSubWindow(v);
    connect(v,SIGNAL(destroyed()),this,SLOT(removeView()));
    s->setSizePolicy(QSizePolicy::Fixed,QSizePolicy::Fixed);
    v->show();
}

void MainWindow::viewObj()
{
    ObjectView* v = new ObjectView(objects_);
    v->setAttribute(Qt::WA_DeleteOnClose,true);
    QMdiSubWindow* s = ui->mdiArea->addSubWindow(v);
    connect(v,SIGNAL(destroyed()),this,SLOT(removeView()));
    s->setSizePolicy(QSizePolicy::Fixed,QSizePolicy::Fixed);
    v->show();
}

void MainWindow::showSVColor()
{
    if(inputs_.empty())return;
    std::vector<MeshBundle<DefaultMesh>::Ptr>::iterator iter;
    for(iter=inputs_.begin();iter!=inputs_.end();++iter)
    {
        if(!*iter||0==(*iter).use_count())return;
        MeshBundle<DefaultMesh>& m = **iter;
        if(m.graph_.voxel_colors.size()==0)return;
        arma::Mat<uint8_t> cmat(
                    (uint8_t*)m.custom_color_.vertex_colors(),
                    4,
                    m.mesh_.n_vertices(),
                    false,
                    true
        );
        m.graph_.sv2pix(m.graph_.voxel_colors,cmat);
    }
}

MainWindow::~MainWindow()
{
    delete ui;
}
