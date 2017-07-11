#include "mainwindow.h"
#include "ui_mainwindow.h"
#include <QFileDialog>
#include <QFileInfo>
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Tools/Utils/Timer.hh>
#include "MeshPairViewerWidget.h"
#include "labspace.h"
#include <vector>
#include <QMdiSubWindow>
#include <QDebug>
#include "objectview.h"
#include "featureview.h"
#include <typeinfo>
#include <fstream>
#include <strstream>
#include "robustcut.h"
#include "iocore.h"
#include "spectrum.h"
#include <QTime>
#include "cube.h"
#include <QDate>
MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow),
    edit_thread_(NULL),edit_widget_(NULL),config_(new Config("./Default.config"))
{
    ui->setupUi(this);

    gl_timer.stop();
    if(!config_->has("Configure"))config_->reload("../Default.config");
    if(!config_->has("Configure"))
    {
        QString msg = "Please Mannually Configure\n";
        QMessageBox::critical(this, windowTitle(), msg);
    }
    connect(ui->actionConfigure,SIGNAL(triggered(bool)),this,SLOT(configure()));
    connect(ui->actionOpen_Inputs,SIGNAL(triggered(bool)),this,SLOT(open_inputs()));
    connect(ui->actionSave_Aligned,SIGNAL(triggered(bool)),this,SLOT(save_aligned()));
    connect(ui->actionSave_Supervoxels,SIGNAL(triggered(bool)),this,SLOT(save_supervoxels()));
    connect(ui->actionLoad_Supervoxels,SIGNAL(triggered(bool)),this,SLOT(load_supervoxels()));
    connect(ui->actionSave_Segments,SIGNAL(triggered(bool)),this,SLOT(save_labels()));
    connect(ui->actionLoad_Segments,SIGNAL(triggered(bool)),this,SLOT(load_labels()));
    connect(ui->actionSave_Base_Segments,SIGNAL(triggered(bool)),this,SLOT(save_base_segs()));
    connect(ui->actionLoad_Base_Segments,SIGNAL(triggered(bool)),this,SLOT(load_base_segs()));
    connect(ui->actionSave_Object_Model,SIGNAL(triggered(bool)),this,SLOT(save_objects()));
    connect(ui->actionLoad_Objects,SIGNAL(triggered(bool)),this,SLOT(load_objects()));
    connect(ui->actionSave_Clusters,SIGNAL(triggered(bool)),this,SLOT(save_cluster()));
    connect(ui->actionLoad_Clusters,SIGNAL(triggered(bool)),this,SLOT(load_cluster()));
    connect(ui->actionSave_Scenes,SIGNAL(triggered(bool)),this,SLOT(save_scenes()));
    connect(ui->actionSave_Points_Index_Picked,SIGNAL(triggered(bool)),this,SLOT(save_pts_index_picked()));
    connect(ui->actionLoad_Points_Index_Picked,SIGNAL(triggered(bool)),this,SLOT(load_pts_index_picked()));
    connect(ui->actionSave_Voxel_Index_Picked,SIGNAL(triggered(bool)),this,SLOT(save_vox_index_picked()));
    connect(ui->actionLoad_Voxel_Index_Picked,SIGNAL(triggered(bool)),this,SLOT(load_vox_index_picked()));
    connect(ui->actionSave_Pix_Order_Functor,SIGNAL(triggered(bool)),this,SLOT(save_Pix_Order_Functor()));
    connect(ui->actionSave_Vox_Order_Functor,SIGNAL(triggered(bool)),this,SLOT(save_Vox_Order_Functor()));
    connect(ui->actionSave_Cube_Color,SIGNAL(triggered(bool)),this,SLOT(save_cube_color()));
    connect(ui->actionLoad_Cube_Color,SIGNAL(triggered(bool)),this,SLOT(load_cube_color()));

    connect(ui->actionGlobal_Align,SIGNAL(triggered(bool)),this,SLOT(start_editing()));
    connect(ui->actionExtract_Background,SIGNAL(triggered(bool)),this,SLOT(start_editing()));
    connect(ui->actionRemove_Zero_Label,SIGNAL(triggered(bool)),this,SLOT(remove_zero_label()));
    connect(ui->actionAnnotator,SIGNAL(triggered(bool)),this,SLOT(start_editing()));
    connect(ui->actionSupervoxel,SIGNAL(triggered(bool)),this,SLOT(start_editing()));
    connect(ui->actionRegionGrow,SIGNAL(triggered(bool)),this,SLOT(start_editing()));
    connect(ui->actionAutomatically,SIGNAL(triggered(bool)),this,SLOT(start_editing()));
    connect(ui->actionMannually,SIGNAL(triggered(bool)),this,SLOT(start_editing()));
    connect(ui->actionV0,SIGNAL(triggered(bool)),this,SLOT(start_editing()));
    connect(ui->actionV1,SIGNAL(triggered(bool)),this,SLOT(start_editing()));
    connect(ui->actionUpdate_Cluster_Center,SIGNAL(triggered(bool)),this,SLOT(start_editing()));
    connect(ui->actionGlobal_Graph_Cut,SIGNAL(triggered(bool)),this,SLOT(start_editing()));
    connect(ui->actionIn_Patch_Graph_Cut,SIGNAL(triggered(bool)),this,SLOT(start_editing()));
    connect(ui->actionIterate,SIGNAL(triggered(bool)),this,SLOT(start_editing()));
    connect(ui->actionGet_Compact_Label,SIGNAL(triggered(bool)),this,SLOT(start_editing()));

    connect(ui->actionJRCS_Old,SIGNAL(triggered(bool)),this,SLOT(start_editing()));
    connect(ui->actionJRCS_Init,SIGNAL(triggered(bool)),this,SLOT(start_editing()));
    connect(ui->actionJRCS_Init_SIHKS,SIGNAL(triggered(bool)),this,SLOT(start_editing()));
    connect(ui->actionDebug_HKS_Clustering,SIGNAL(triggered(bool)),this,SLOT(start_editing()));
    connect(ui->actionJRCS_Init_Bernoulli,SIGNAL(triggered(bool)),this,SLOT(start_editing()));
    connect(ui->actionMore_Box,SIGNAL(triggered(bool)),this,SLOT(start_editing()));
    connect(ui->actionJRCS_Opt_Basic,SIGNAL(triggered(bool)),this,SLOT(start_editing()));
    connect(ui->actionJRCS_Opt_AONI,SIGNAL(triggered(bool)),this,SLOT(start_editing()));
    connect(ui->actionJRCS_Opt_AOPT,SIGNAL(triggered(bool)),this,SLOT(start_editing()));
    connect(ui->actionJRCS_Opt_Spectrum,SIGNAL(triggered(bool)),this,SLOT(start_editing()));
    connect(ui->actionJRCS_Opt_Bilateral,SIGNAL(triggered(bool)),this,SLOT(start_editing()));
    connect(ui->actionJRCS_Opt_Primitive,SIGNAL(triggered(bool)),this,SLOT(start_editing()));
    connect(ui->actionJRCS_Opt_Cube,SIGNAL(triggered(bool)),this,SLOT(start_editing()));
    connect(ui->actionJRCS_Opt_Box,SIGNAL(triggered(bool)),this,SLOT(start_editing()));
    connect(ui->actionGo_Over,SIGNAL(triggered(bool)),this,SLOT(goOver()));

    connect(ui->actionGDCoord,SIGNAL(triggered(bool)),this,SLOT(start_editing()));
    connect(ui->actionScene_Maker,SIGNAL(triggered(bool)),this,SLOT(make_scene()));
    connect(ui->actionEstimate_IOU,SIGNAL(triggered(bool)),this,SLOT(calculate_iou()));
    connect(ui->actionEstimate_Registration_for_JRCS,SIGNAL(triggered(bool)),this,SLOT(calculate_fit()));

    connect(ui->actionRegionGrowRGB,SIGNAL(triggered(bool)),this,SLOT(start_editing()));
    connect(ui->actionSort_AGD,SIGNAL(triggered(bool)),this,SLOT(start_editing()));
    connect(ui->actionDebug_Convexity,SIGNAL(triggered(bool)),this,SLOT(start_editing()));
    connect(ui->actionDebug_Color,SIGNAL(triggered(bool)),this,SLOT(start_editing()));
    connect(ui->actionDebug_Dist,SIGNAL(triggered(bool)),this,SLOT(start_editing()));
    connect(ui->actionDebug_W,SIGNAL(triggered(bool)),this,SLOT(start_editing()));
    connect(ui->actionCut_Graph,SIGNAL(triggered(bool)),this,SLOT(start_editing()));
    connect(ui->actionBase_Segments,SIGNAL(triggered(bool)),this,SLOT(start_editing()));
    connect(ui->actionConsensus_Segment,SIGNAL(triggered(bool)),this,SLOT(start_editing()));
    connect(ui->actionShow_Base_Segments,SIGNAL(triggered(bool)),this,SLOT(start_editing()));
    connect(ui->actionMakeIncomplete,SIGNAL(triggered(bool)),this,SLOT(start_editing()));

    connect(ui->actionLab_Color_Space,SIGNAL(triggered(bool)),this,SLOT(showLab()));
    connect(ui->actionObject_View,SIGNAL(triggered(bool)),this,SLOT(viewObj()));
    connect(ui->actionSupervoxel_Color,SIGNAL(triggered(bool)),this,SLOT(showSVColor()));
    connect(ui->actionFeature_View,SIGNAL(triggered(bool)),this,SLOT(showFeature()));
    connect(ui->actionIndex_By_Color,SIGNAL(triggered(bool)),this,SLOT(showIndex()));
    connect(ui->actionSpectral_Function,SIGNAL(triggered(bool)),this,SLOT(showSpectralFunc()));
    connect(ui->actionColor_By_Cube,SIGNAL(triggered(bool)),this,SLOT(custom_color_from_cube()));

    connect(ui->actionLAPACKE_dggsvd,SIGNAL(triggered(bool)),this,SLOT(LAPACKE_dggsvd_test()));
    connect(ui->actionInside_Bounding_Box,SIGNAL(triggered(bool)),this,SLOT(Inside_BBox_test()));
    connect(ui->actionAGD_test,SIGNAL(triggered(bool)),this,SLOT(agd_test()));
    connect(ui->actionJRCS_Plate,SIGNAL(triggered(bool)),this,SLOT(jrcs_plate_test()));
    connect(ui->actionJRCS_Cube,SIGNAL(triggered(bool)),this,SLOT(jrcs_cube_test()));

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

void MainWindow::open_inputs(QDir& dir)
{
    QStringList fileNames = dir.entryList();
    QStringList filepath;
    foreach(QString name,fileNames)
    {
        filepath << dir.absoluteFilePath(name);
    }
    if (!filepath.isEmpty())
    {
        open_inputs(filepath);
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
    if(config_->has("ascii"))
    {
        opt.unset(OpenMesh::IO::Options::Binary);
    }
    std::string path = dirName.toStdString();
    MeshBundle<DefaultMesh>::PtrList::iterator iter;
    for(iter=inputs_.begin();iter!=inputs_.end();++iter)
    {
        MeshBundle<DefaultMesh>::Ptr ptr = *iter;
        if(!OpenMesh::IO::write_mesh(ptr->mesh_,path+"/"+ptr->name_+".ply",opt,13)){
            std::cerr<<"can't save to:"<<path+"/"+ptr->name_+".ply"<<std::endl;
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

void MainWindow::view_input(QWidget* widget)
{
    gl_timer.stop();
    connect(ui->actionCustom_Color,SIGNAL(toggled(bool)),widget,SLOT(use_custom_color(bool)));
    connect(&gl_timer,SIGNAL(timeout()),widget,SLOT(updateGL()));
    mesh_views_.push_back(WidgetPtr(widget));
    showInMdi((QWidget*)widget,Qt::Widget|Qt::WindowMinMaxButtonsHint);
    gl_timer.start(100);
}

void MainWindow::save_labels(QString dirName)
{
    if(dirName.isEmpty())
    {
        dirName = QFileDialog::getExistingDirectory(
                this,
                tr("Save Labels"),
                tr("../Dev_Data/")
                );
    }
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
        if(config_->has("ascii"))
        {
            if(!iter->save(filepath.toStdString(),arma::raw_ascii))
            {
                QString msg = "Failed to Save "+filepath+"\n";
                QMessageBox::critical(this, windowTitle(), msg);
                return;
            }
        }else{
            if(!iter->save(filepath.toStdString(),arma::arma_binary))
            {
                QString msg = "Failed to Save "+filepath+"\n";
                QMessageBox::critical(this, windowTitle(), msg);
                return;
            }
        }
        if(miter==inputs_.end())break;
        ++miter;
    }
}

void MainWindow::save_cube_color(QString fileName)
{
    if(fileName.isEmpty())
    {
        fileName = QFileDialog::getSaveFileName(this,
                                                tr("Save Cube Color"),
                                                tr("./"),
                                                tr("cube color(*.arma)"));
    }
    if (fileName.isEmpty())return;
    arma::uvec cube_color(Common::Cube::color_size());
    for(int i=0;i<cube_color.size();++i)
    {
        cube_color(i) = Common::Cube::colorFromLabel(i+1);
    }
    if(config_->has("ascii"))
    {
        if(!cube_color.save(fileName.toStdString(),arma::raw_ascii))
        {
            QString msg = "Failed to Save "+fileName+"\n";
            QMessageBox::critical(this, windowTitle(), msg);
            return;
        }
    }else{
        if(!cube_color.save(fileName.toStdString(),arma::arma_binary))
        {
            QString msg = "Failed to Save "+fileName+"\n";
            QMessageBox::critical(this, windowTitle(), msg);
            return;
        }
    }

}

void MainWindow::load_cube_color(QString fileName)
{
    if(fileName.isEmpty())
    {
        fileName = QFileDialog::getOpenFileName(this,
                                                tr("Load Cube Color"),
                                                tr("./"),
                                                tr("cube color(*.arma)"));
    }
    if (fileName.isEmpty())return;
    arma::uvec cube_color;
    if(!cube_color.load(fileName.toStdString()))
    {
        QString msg = "Failed to Load "+fileName+"\n";
        QMessageBox::critical(this, windowTitle(), msg);
        return;
    }
    Common::Cube::set_color(cube_color);
    if(labels_.size()!=inputs_.size())
    {
        return;
    }
    std::vector<arma::uvec>::iterator iter;
    std::vector<MeshBundle<DefaultMesh>::Ptr>::iterator miter;
    miter = inputs_.begin();
    for(iter=labels_.begin();iter!=labels_.end();++iter)
    {
        MeshBundle<DefaultMesh>& m = **miter;
        Common::Cube::colorByLabel((uint32_t*)m.custom_color_.vertex_colors(),m.custom_color_.size(),*iter);
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

void MainWindow::custom_color_from_cube()
{
    if(labels_.size()!=inputs_.size())
    {
        return;
    }
    std::vector<arma::uvec>::iterator iter;
    std::vector<MeshBundle<DefaultMesh>::Ptr>::iterator miter;
    miter = inputs_.begin();
    Common::Cube::reset_color_set();
    for(iter=labels_.begin();iter!=labels_.end();++iter)
    {
        MeshBundle<DefaultMesh>& m = **miter;
        Common::Cube::colorByLabel((uint32_t*)m.custom_color_.vertex_colors(),m.custom_color_.size(),*iter);
        if(miter==inputs_.end())break;
        ++miter;
    }
}

void MainWindow::save_base_segs(QString dirName)
{
    if(dirName.isEmpty())
    {
        dirName = QFileDialog::getExistingDirectory(
                this,
                tr("Save Base Segments"),
                tr("../Dev_Data/")
                );
    }
    if(dirName.isEmpty())return;
    std::vector<arma::umat>::iterator iter;
    std::vector<MeshBundle<DefaultMesh>::Ptr>::iterator miter;
    QDir dir;
    dir.setPath(dirName);
    miter = inputs_.begin();
    for(iter=RobustCut::base_segment_list_.begin();iter!=RobustCut::base_segment_list_.end();++iter)
    {
        QString filepath = dir.absoluteFilePath(
                    QString::fromStdString((*miter)->name_+".umat.arma")
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

void MainWindow::load_base_segs()
{
    if(inputs_.empty())
    {
        QString msg = "Please Load Inputs First\n";
        QMessageBox::critical(this, windowTitle(), msg);
    }
    QString dirName = QFileDialog::getExistingDirectory(
                this,
                tr("Load Base Segments"),
                tr("../Dev_Data/")
                );
    if(dirName.isEmpty())return;
    std::vector<arma::umat>::iterator iter;
    std::vector<MeshBundle<DefaultMesh>::Ptr>::iterator miter;
    QDir dir;
    dir.setPath(dirName);
    miter = inputs_.begin();
    if(RobustCut::base_segment_list_.size()!=inputs_.size())
    {
        RobustCut::base_segment_list_.resize(inputs_.size());
    }
    for(iter=RobustCut::base_segment_list_.begin();iter!=RobustCut::base_segment_list_.end();++iter)
    {
        arma::umat label;
        QString filepath = dir.absoluteFilePath(
                    QString::fromStdString((*miter)->name_+".umat.arma")
                    );
        if(!label.load(filepath.toStdString()))
        {
            QString msg = "Failed to Load "+filepath+"\n";
            QMessageBox::critical(this, windowTitle(), msg);
            return;
        }
        *iter = label;
        if(miter==inputs_.end())break;
        ++miter;
    }
}

void MainWindow::save_objects(QString dirName)
{
    if(dirName.isEmpty())
    {
        dirName = QFileDialog::getExistingDirectory(
                this,
                tr("Save Objects"),
                tr("../Dev_Data/")
                );
    }
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

void MainWindow::save_supervoxels(QString dirName)
{
    if(dirName.isEmpty())
    {
        dirName = QFileDialog::getExistingDirectory(
                this,
                tr("Save Supervoxels"),
                tr("../Dev_Data/")
                );
    }
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

void MainWindow::save_cluster(QString dirName)
{
    std::cerr<<"saving cluster"<<std::endl;
    if(dirName.isEmpty())
    {
        dirName = QFileDialog::getExistingDirectory(
                this,
                tr("Save Cluster"),
                tr("../Dev_Data/")
                );
    }
    if(dirName.isEmpty())return;
    QDir dir;
    dir.setPath(dirName);
    QString filepath = dir.absoluteFilePath(
                tr("Base.fmat.arma")
                );
    if(!feature_base_.save(filepath.toStdString(),arma::arma_binary))
    {
        QString msg = "Failed to Save "+filepath+"\n";
        QMessageBox::critical(this, windowTitle(), msg);
        return;
    }
    filepath = dir.absoluteFilePath(
                tr("Center.fmat.arma")
                );
    if(!feature_centers_.save(filepath.toStdString(),arma::arma_binary))
    {
        QString msg = "Failed to Save "+filepath+"\n";
        QMessageBox::critical(this, windowTitle(), msg);
        return;
    }
}

void MainWindow::load_cluster()
{
    QString dirName = QFileDialog::getExistingDirectory(
            this,
            tr("Load Cluster"),
            tr("../Dev_Data/")
            );
    if(dirName.isEmpty())return;
    QDir dir;
    dir.setPath(dirName);
    QString filepath = dir.absoluteFilePath(
                QString::fromStdString("Base.fmat.arma")
                );
    if(!feature_base_.load(filepath.toStdString()))
    {
        QString msg = "Failed to Load "+filepath+"\n";
        QMessageBox::critical(this, windowTitle(), msg);
        return;
    }
    filepath = dir.absoluteFilePath(
                QString::fromStdString("Center.fmat.arma")
                );
    if(!feature_centers_.load(filepath.toStdString()))
    {
        QString msg = "Failed to Load "+filepath+"\n";
        QMessageBox::critical(this, windowTitle(), msg);
        return;
    }
}

void MainWindow::save_vox_index_picked(QString dirName)
{
    if(inputs_.empty()){
        QString msg = "Load inputs first";
        QMessageBox::warning(this, windowTitle(), msg);
        return;
    }
    if(inputs_[0]->graph_.empty())
    {
        QString msg = "Load supervoxels first";
        QMessageBox::warning(this, windowTitle(), msg);
        return;
    }
    if(dirName.isEmpty())dirName = QFileDialog::getExistingDirectory(
            this,
            tr("Save Index of Picked Points"),
            tr("./debug/HKS/")
            );
    if(dirName.isEmpty())return;
    std::vector<WidgetPtr>::iterator iter;
    arma::uword index = 1;
    for(iter=mesh_views_.begin();iter!=mesh_views_.end();++iter)
    {
        MeshPairViewerWidget* w = qobject_cast<MeshPairViewerWidget*>(*iter);
        if(w)
        {
            if( index > inputs_.size())
            {
                std::cerr<<"the number of input doesn't match the view number"<<std::endl;
                return;
            }
            MeshBundle<DefaultMesh>& mesh = *inputs_[index-1];
            QString path;
            path = path.sprintf("vox_index%02u.mat",index);
            if(!w->first_selected().empty())
            {
                arma::uvec selected;
                mesh.graph_.getSvIndex(w->first_selected(),selected);
                MATIO::save_to_matlab(selected,(dirName+"/"+path).toStdString(),"X");
            }
            ++index;
        }
    }
}

void MainWindow::load_vox_index_picked()
{
    if(inputs_.empty()){
        QString msg = "Load inputs first";
        QMessageBox::warning(this, windowTitle(), msg);
        return;
    }
    if(inputs_[0]->graph_.empty())
    {
        QString msg = "Load supervoxels first";
        QMessageBox::warning(this, windowTitle(), msg);
        return;
    }
    QString dirName = QFileDialog::getExistingDirectory(
            this,
            tr("Load Index of Picked Voxels"),
            tr("./debug/HKS/")
            );
    if(dirName.isEmpty())return;
    std::vector<WidgetPtr>::iterator iter;
    arma::uword index = 1;
    QDir dir;
    dir.setPath(dirName);
    for(iter=mesh_views_.begin();iter!=mesh_views_.end();++iter)
    {
        MeshPairViewerWidget* w = qobject_cast<MeshPairViewerWidget*>(*iter);
        if(w)
        {
            if( index > inputs_.size())
            {
                std::cerr<<"the number of input doesn't match the view number"<<std::endl;
                return;
            }
            MeshBundle<DefaultMesh>& mesh = *inputs_[index-1];
            disconnect(&gl_timer,SIGNAL(timeout()),w,SLOT(updateGL()));
            QString path;
            path = path.sprintf("vox_index%02u.mat",index);
            QFileInfo info(dir.absoluteFilePath(path));
            if(info.exists())
            {
                arma::uvec selected;
                MATIO::load_to_arma(selected,info.absoluteFilePath().toStdString(),"X");
                arma::uvec pix;
                mesh.graph_.getPixIndex(selected,pix);
                w->first_selected().resize(pix.size());
                w->first_selected() = arma::conv_to<std::vector<arma::uword>>::from(pix);
            }else{
                std::cerr<<"Failed to find:"<<info.absoluteFilePath().toStdString()<<std::endl;
            }
            ++index;
            connect(&gl_timer,SIGNAL(timeout()),w,SLOT(updateGL()));
        }
    }

}

void MainWindow::save_pts_index_picked(QString dirName)
{
    if(dirName.isEmpty())dirName = QFileDialog::getExistingDirectory(
            this,
            tr("Save Index of Picked Points"),
            tr("./debug/HKS/")
            );
    if(dirName.isEmpty())return;
    std::vector<WidgetPtr>::iterator iter;
    arma::uword index = 1;
    for(iter=mesh_views_.begin();iter!=mesh_views_.end();++iter)
    {
        MeshPairViewerWidget* w = qobject_cast<MeshPairViewerWidget*>(*iter);
        if(w)
        {
            QString path;
            path = path.sprintf("pts_index%02u.mat",index);
            if(!w->first_selected().empty())
            {
                arma::uvec selected(w->first_selected());
                MATIO::save_to_matlab(selected,(dirName+"/"+path).toStdString(),"X");
            }
            ++index;
        }
    }
}

void MainWindow::save_Pix_Order_Functor(QString dirName)
{
    if(dirName.isEmpty())dirName = QFileDialog::getExistingDirectory(
            this,
            tr("Save Pix Order Functor"),
            tr("../Dev_Data/")
            );
    if(dirName.isEmpty())return;
    std::vector<WidgetPtr>::iterator iter;
    arma::uword index = 1;
    for(iter=mesh_views_.begin();iter!=mesh_views_.end();++iter)
    {
        MeshPairViewerWidget* w = qobject_cast<MeshPairViewerWidget*>(*iter);
        if(w)
        {
            QString path;
            path = path.sprintf("FuncOn%sPix_",w->first_ptr()->name_.c_str());
            path += QDate::currentDate().toString(tr("yyyyMMdd"))+QTime::currentTime().toString(tr("hhmmsszzz"));
            path += tr(".mat");
            if(!w->first_selected().empty())
            {
                arma::uvec selected(w->first_selected());
                arma::vec value(w->first_ptr()->mesh_.n_vertices());
                value.fill(selected.size()+1);
                value(selected) = arma::linspace<arma::vec>(1,selected.size(),selected.size());
                MATIO::save_to_matlab(value,(dirName+"/"+path).toStdString(),"X");
            }
            ++index;
        }
    }
}

void MainWindow::save_Vox_Order_Functor(QString dirName)
{
    if(inputs_.empty()){
        QString msg = "Load inputs first";
        QMessageBox::warning(this, windowTitle(), msg);
        return;
    }
    if(inputs_[0]->graph_.empty())
    {
        QString msg = "Load supervoxels first";
        QMessageBox::warning(this, windowTitle(), msg);
        return;
    }
    if(dirName.isEmpty())dirName = QFileDialog::getExistingDirectory(
            this,
            tr("Save Vox Order Functor"),
            tr("../Dev_Data/")
            );
    if(dirName.isEmpty())return;
    std::vector<WidgetPtr>::iterator iter;
    arma::uword index = 1;
    for(iter=mesh_views_.begin();iter!=mesh_views_.end();++iter)
    {
        MeshPairViewerWidget* w = qobject_cast<MeshPairViewerWidget*>(*iter);
        if(w)
        {
            QString path;
            path = path.sprintf("FuncOn%sVox_",w->first_ptr()->name_.c_str());
            path += QDate::currentDate().toString(tr("yyyyMMdd"))+QTime::currentTime().toString(tr("hhmmsszzz"));
            path += tr(".mat");
            if(!w->first_selected().empty())
            {
                arma::uvec selected(w->first_selected());
                arma::uvec svSelected;
                w->first_ptr()->graph_.getSvIndex(selected,svSelected);
                arma::vec vox_value(w->first_ptr()->graph_.size(),arma::fill::zeros);
                svSelected -= 1;
                vox_value(svSelected) = arma::linspace<arma::vec>(1,svSelected.size(),svSelected.size());
                arma::vec value(w->first_ptr()->mesh_.n_vertices(),arma::fill::zeros);
                arma::uvec indices = w->first_ptr()->graph_.voxel_label;
                arma::uvec idx = arma::find( indices > 0 );
                indices(idx) -= 1;
                value = vox_value(indices);
                MATIO::save_to_matlab(value,(dirName+"/"+path).toStdString(),"X");
            }
            ++index;
        }
    }
}

void MainWindow::load_pts_index_picked()
{
    QString dirName = QFileDialog::getExistingDirectory(
            this,
            tr("Load Index of Picked Points"),
            tr("./debug/HKS/")
            );
    if(dirName.isEmpty())return;
    std::vector<WidgetPtr>::iterator iter;
    arma::uword index = 1;
    QDir dir;
    dir.setPath(dirName);
    for(iter=mesh_views_.begin();iter!=mesh_views_.end();++iter)
    {
        MeshPairViewerWidget* w = qobject_cast<MeshPairViewerWidget*>(*iter);
        if(w)
        {
            disconnect(&gl_timer,SIGNAL(timeout()),w,SLOT(updateGL()));
            QString path;
            path = path.sprintf("pts_index%02u.mat",index);
            QFileInfo info(dir.absoluteFilePath(path));
            if(info.exists())
            {
                arma::uvec selected;
                MATIO::load_to_arma(selected,info.absoluteFilePath().toStdString(),"X");
//                std::cerr<<selected<<std::endl;
                w->first_selected().resize(selected.size());
                w->first_selected() = arma::conv_to<std::vector<arma::uword>>::from(selected);
            }else{
                std::cerr<<"Failed to find:"<<info.absoluteFilePath().toStdString()<<std::endl;
            }
            ++index;
            connect(&gl_timer,SIGNAL(timeout()),w,SLOT(updateGL()));
        }
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
        stream<<path<<"/layouts/obj"<<index<<".obj";
        stream>>objfilename;
        objfile.open(objfilename,objfile.out);
        if(!objfile.is_open())
        {
            std::cerr<<"Failed to Open:"<<objfilename<<std::endl;
        }
        objfile<<"#obj obj"<<index<<std::endl;
        arma::fmat layout_mat;
        if(obj_ptr->layout(layout_mat,-1));
        {
            layout_stream.clear();
            layout_stream.str("");
            size_t cnt = 1;
            for(size_t c=0;c<layout_mat.n_cols;++c)
            {
                layout_stream<<"v "<<layout_mat(0,c)
                            <<" "<<layout_mat(1,c)
                           <<" "<<layout_mat(2,c)
                          <<std::endl;
            }
            objfile << layout_stream.str();
            objfile<<"f "<<cnt+0<<" "<<cnt+3<<" "<<cnt+2<<" "<<cnt+1<<std::endl;
            objfile<<"f "<<cnt+4<<" "<<cnt+5<<" "<<cnt+6<<" "<<cnt+7<<std::endl;
            objfile<<"f "<<cnt+0<<" "<<cnt+4<<" "<<cnt+7<<" "<<cnt+3<<std::endl;
            objfile<<"f "<<cnt+0<<" "<<cnt+1<<" "<<cnt+5<<" "<<cnt+4<<std::endl;
            objfile<<"f "<<cnt+1<<" "<<cnt+2<<" "<<cnt+6<<" "<<cnt+5<<std::endl;
            objfile<<"f "<<cnt+2<<" "<<cnt+3<<" "<<cnt+7<<" "<<cnt+6<<std::endl;
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
    OpenMesh::IO::Options opt;
    opt+=OpenMesh::IO::Options::Binary;
    opt+=OpenMesh::IO::Options::VertexColor;
    for(size_t fIndex = 0 ; fIndex < inputs_.size() ; ++fIndex )
    {
        stream.clear();
        stream.str("");
        stream<<path<<"/layouts/scene_"<<fIndex;
        stream>>objfilename;
        objfile.open(objfilename+".obj",objfile.out);
        DefaultMesh scene_layout;
        scene_layout.request_vertex_colors();
        MeshColor<DefaultMesh> scene_layout_color(scene_layout);
        std::vector<arma::uword> labels;
        std::vector<DefaultMesh::VertexHandle> vhandles;
        if(!objfile.is_open())
        {
            std::cerr<<"Failed to Open:"<<objfilename+".obj"<<std::endl;
        }

        std::vector<ObjModel::Ptr>::iterator iter;
        //counting objects number
        size_t obj_number = 0;
        for(iter=objects_.begin();iter!=objects_.end();++iter)
        {
            arma::fmat layout_mat;
            if( ( (*iter)->layout(layout_mat,fIndex) ) )
            {
                ++obj_number;
            }
        }
        objfile<<"#obj_number "<<obj_number<<std::endl;
        size_t cnt = 1;
        size_t oIndex = 1;
        for(iter=objects_.begin();iter!=objects_.end();++iter)
        {
            arma::fmat layout_mat;
            if( ( (*iter)->layout(layout_mat,fIndex) ) )
            {
                objfile<<"#obj obj"<<oIndex<<std::endl;
                for(size_t c=0;c<layout_mat.n_cols;++c)
                {
                    objfile<<"v "<<layout_mat(0,c)
                          <<" "<<layout_mat(1,c)
                         <<" "<<layout_mat(2,c)
                        <<std::endl;
                    vhandles.push_back(
                                scene_layout.add_vertex(
                                    DefaultMesh::Point(
                                        layout_mat(0,c),
                                        layout_mat(1,c),
                                        layout_mat(2,c)
                                        )
                                    )
                    );
                    labels.push_back(oIndex);
                }
                objfile<<"f "<<cnt+0<<" "<<cnt+3<<" "<<cnt+2<<" "<<cnt+1<<std::endl;
                objfile<<"f "<<cnt+4<<" "<<cnt+5<<" "<<cnt+6<<" "<<cnt+7<<std::endl;
                objfile<<"f "<<cnt+0<<" "<<cnt+4<<" "<<cnt+7<<" "<<cnt+3<<std::endl;
                objfile<<"f "<<cnt+0<<" "<<cnt+1<<" "<<cnt+5<<" "<<cnt+4<<std::endl;
                objfile<<"f "<<cnt+1<<" "<<cnt+2<<" "<<cnt+6<<" "<<cnt+5<<std::endl;
                objfile<<"f "<<cnt+2<<" "<<cnt+3<<" "<<cnt+7<<" "<<cnt+6<<std::endl;
                {
                    vhandles.clear();
                    vhandles.push_back(scene_layout.vertex_handle(cnt-1+0));
                    vhandles.push_back(scene_layout.vertex_handle(cnt-1+3));
                    vhandles.push_back(scene_layout.vertex_handle(cnt-1+2));
                    vhandles.push_back(scene_layout.vertex_handle(cnt-1+1));
                    scene_layout.add_face(vhandles);

                    vhandles.clear();
                    vhandles.push_back(scene_layout.vertex_handle(cnt-1+4));
                    vhandles.push_back(scene_layout.vertex_handle(cnt-1+5));
                    vhandles.push_back(scene_layout.vertex_handle(cnt-1+6));
                    vhandles.push_back(scene_layout.vertex_handle(cnt-1+7));
                    scene_layout.add_face(vhandles);

                    vhandles.clear();
                    vhandles.push_back(scene_layout.vertex_handle(cnt-1+0));
                    vhandles.push_back(scene_layout.vertex_handle(cnt-1+4));
                    vhandles.push_back(scene_layout.vertex_handle(cnt-1+7));
                    vhandles.push_back(scene_layout.vertex_handle(cnt-1+3));
                    scene_layout.add_face(vhandles);

                    vhandles.clear();
                    vhandles.push_back(scene_layout.vertex_handle(cnt-1+0));
                    vhandles.push_back(scene_layout.vertex_handle(cnt-1+1));
                    vhandles.push_back(scene_layout.vertex_handle(cnt-1+5));
                    vhandles.push_back(scene_layout.vertex_handle(cnt-1+4));
                    scene_layout.add_face(vhandles);

                    vhandles.clear();
                    vhandles.push_back(scene_layout.vertex_handle(cnt-1+1));
                    vhandles.push_back(scene_layout.vertex_handle(cnt-1+2));
                    vhandles.push_back(scene_layout.vertex_handle(cnt-1+6));
                    vhandles.push_back(scene_layout.vertex_handle(cnt-1+5));
                    scene_layout.add_face(vhandles);

                    vhandles.clear();
                    vhandles.push_back(scene_layout.vertex_handle(cnt-1+2));
                    vhandles.push_back(scene_layout.vertex_handle(cnt-1+3));
                    vhandles.push_back(scene_layout.vertex_handle(cnt-1+7));
                    vhandles.push_back(scene_layout.vertex_handle(cnt-1+6));
                    scene_layout.add_face(vhandles);
                }
                cnt+=8;
            }
            ++oIndex;
        }
        arma::uvec scene_label(labels);
        scene_layout_color.fromlabel(scene_label);
        arma::Mat<uint8_t> smat((uint8_t*)scene_layout.vertex_colors(),3,scene_layout.n_vertices(),false,true);
        arma::Mat<uint8_t> cmat((uint8_t*)scene_layout_color.vertex_colors(),4,scene_layout.n_vertices(),false,true);
        smat = cmat.rows(0,2);
        objfile.close();
        OpenMesh::IO::write_mesh(scene_layout,objfilename+".ply",opt,13);
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

void MainWindow::showSpectralFunc()
{
    if(inputs_.empty()){
        ui->statusBar->showMessage(QString("! No Input"));
        return;
    }
    if(inputs_[0]->graph_.empty()){
        ui->statusBar->showMessage(QString("! No Supervoxel"));
        return;
    }
    Spectrum* v = new Spectrum(inputs_,config_);
    v->setAttribute(Qt::WA_DeleteOnClose,true);
    QMdiSubWindow* s = ui->mdiArea->addSubWindow(v);
    connect(v,SIGNAL(destroyed()),this,SLOT(removeView()));
    s->setSizePolicy(QSizePolicy::Minimum,QSizePolicy::Minimum);
    v->show();
}

void MainWindow::viewObj()
{
    if(objects_.empty()){
        ui->statusBar->showMessage(QString("! No Objects"));
        return;
    }
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

void MainWindow::showFeature()
{
    if(inputs_.empty())return;
    featureview* v = new featureview(inputs_,labels_,feature_base_,feature_centers_);
    if(!v->configure(config_))return;
    v->setAttribute(Qt::WA_DeleteOnClose,true);
    QMdiSubWindow* s = ui->mdiArea->addSubWindow(v);
    connect(v,SIGNAL(destroyed()),this,SLOT(removeView()));
    v->init();
    v->show();
}

void MainWindow::showIndex()
{
    if(inputs_.empty())return;
    std::vector<MeshBundle<DefaultMesh>::Ptr>::iterator iter;
    for(iter=inputs_.begin();iter!=inputs_.end();++iter)
    {
        if(!*iter||0==(*iter).use_count())continue;
        MeshBundle<DefaultMesh>& m = **iter;
        m.custom_color_.fromIndex();
    }
    ui->actionCustom_Color->setChecked(true);
}

void MainWindow::calculate_iou(QString dir0,QString dir1)
{
    if(inputs_.empty())
    {
        QString msg = "Please Load Inputs First\n";
        QMessageBox::critical(this, windowTitle(), msg);
    }
    if(dir0.isEmpty())
    {
        dir0 = QFileDialog::getExistingDirectory(
                this,
                tr("Load Labels 0"),
                tr("../Dev_Data/")
                );
    }
    if(dir0.isEmpty())return;

    if(dir1.isEmpty())
    {
        dir1 = QFileDialog::getExistingDirectory(
                this,
                tr("Load Labels 1"),
                tr("../Dev_Data/")
                );
    }
    if(dir1.isEmpty())return;

    QDir dir;
    QString filepath;
    float iou = 0.0;
    float min = std::numeric_limits<float>::max();
    float max = std::numeric_limits<float>::lowest();
    int mini,maxi;
    for(int i=0;i<inputs_.size();++i)
    {
        dir.setPath(dir0);
        arma::uvec lbl0,lbl1;
        filepath = dir.absoluteFilePath(
                    QString::fromStdString((inputs_[i])->name_+".label.arma")
                    );
        if(!lbl0.load(filepath.toStdString()))
        {
            QString msg = "Failed to Load "+filepath+"\n";
            QMessageBox::critical(this, windowTitle(), msg);
            return;
        }
        if( lbl0.size() != (inputs_[i])->mesh_.n_vertices() )
        {
            QString msg = "Some size of this set of labels doesn't match the inputs\n";
            QMessageBox::critical(this, windowTitle(), msg);
            return;
        }
        dir.setPath(dir1);
        filepath = dir.absoluteFilePath(
                    QString::fromStdString((inputs_[i])->name_+".label.arma")
                    );
        if(!lbl1.load(filepath.toStdString()))
        {
            QString msg = "Failed to Load "+filepath+"\n";
            QMessageBox::critical(this, windowTitle(), msg);
            return;
        }
        if( lbl1.size() != (inputs_[i])->mesh_.n_vertices() )
        {
            QString msg = "Some size of this set of labels doesn't match the inputs\n";
            QMessageBox::critical(this, windowTitle(), msg);
            return;
        }
        arma::uvec intersect = arma::find(lbl0==lbl1);
        iou = float(intersect.size());
        iou = iou / ( 2.0*float(lbl0.size()) - iou );
        std::cout<<(inputs_[i])->name_<<":"<<iou<<std::endl;
        if(min>iou)
        {
            mini = i;
            min=iou;
        }
        if(max<iou)
        {
            maxi = i;
            max=iou;
        }
    }
    std::cout<<"min:"<<(inputs_[mini])->name_<<":"<<min<<std::endl;
    std::cout<<"max:"<<(inputs_[maxi])->name_<<":"<<max<<std::endl;
    std::cout.flush();
}

void MainWindow::calculate_fit(QString dirName)
{
    if(inputs_.empty())
    {
        QString msg = "Please Load Inputs First\n";
        QMessageBox::critical(this, windowTitle(), msg);
    }
    if(dirName.isEmpty())
    {
        dirName = QFileDialog::getExistingDirectory(
                this,
                tr("Load RT"),
                tr("../Dev_Data/")
                );
    }
    if(dirName.isEmpty())return;
    std::vector<arma::uvec>::iterator iter;
    std::vector<MeshBundle<DefaultMesh>::Ptr>::iterator miter = inputs_.begin();
    QDir dir;
    dir.setPath(dirName);
    arma::fmat v0;
    arma::fmat v1((float*)(inputs_.back()->mesh_.points()),3,inputs_.back()->mesh_.n_vertices(),true,true);
    float min = std::numeric_limits<float>::max();
    float max = std::numeric_limits<float>::lowest();
    int mini,maxi;
    int i=0;
    for(iter=labels_.begin();iter!=labels_.end();++iter)
    {
        QString filepath = dir.absoluteFilePath(
                    QString::fromStdString((*miter)->name_+".txt")
                    );
        std::ifstream in;
        in.open(filepath.toStdString());
        v0 = v1;
        v1 = arma::fmat((float*)(*miter)->mesh_.points(),3,(*miter)->mesh_.n_vertices(),true,true);
        arma::fmat R(3,3);
        arma::fvec t(3);
        arma::uvec & label = *iter;
        arma::uword maxl = arma::max(label);
        if( maxl==0 || label.empty() )
        {
            in >> R(0,0) >> R(0,1) >> R(0,2) >> t(0);
            in >> R(1,0) >> R(1,1) >> R(1,2) >> t(1);
            in >> R(2,0) >> R(2,1) >> R(2,2) >> t(2);
            std::cerr<<i<<std::endl;
            std::cerr<<R<<std::endl;
            v0 = R*v0;
            v0.each_col() += t;
        }else{
            for(arma::uword l=arma::min(label);l<maxl;++l)
            {
                in >> R(0,0) >> R(0,1) >> R(0,2) >> t(0);
                in >> R(1,0) >> R(1,1) >> R(1,2) >> t(1);
                in >> R(2,0) >> R(2,1) >> R(2,2) >> t(2);
                arma::uvec indices = arma::find(label==l);
                arma::fmat cols = v0.cols(indices);
                cols = R*cols;
                cols.each_col() += t;
                v0.cols(indices) = cols;
            }
        }
        in.close();
        arma::fvec scale = ( arma::max(v1,1) - arma::min(v1,1) );
        float s = arma::mean(scale);
        arma::fmat e = arma::square(v0 - v1);
        arma::frowvec e2 = arma::sqrt(arma::sum(e));
        float error = arma::mean( e2 );
        error /= s;
        std::cout<<(*miter)->name_<<":"<<error<<std::endl;
        if(min>error)
        {
            mini = i;
            min=error;
        }
        if(max<error)
        {
            maxi = i;
            max=error;
        }
        if(miter==inputs_.end())break;
        ++miter;
        ++i;
    }
    std::cout<<"min:"<<(inputs_[mini])->name_<<":"<<min<<std::endl;
    std::cout<<"max:"<<(inputs_[maxi])->name_<<":"<<max<<std::endl;
    std::cout.flush();
}

void MainWindow::goOver()
{
    QString fileName = QFileDialog::getOpenFileName(this,
        tr("Go Over"),
        tr("./"),
        tr("GoOver (*.goover);;"
        "All Files (*)"));
    if(fileName.isEmpty())return;
    QFile inFile(fileName);
    inFile.open(inFile.ReadOnly);
    QTextStream stream(&inFile);
    QTextStream configtxt;
    Config::Ptr config;
    config.reset(new Config());
    uint32_t num = 0;
    stream >> num;
    QString txt;
    stream >> txt;
    QAction* act = getActionByText(txt);
    if(!act){
        std::cerr<<"Can't find the action in this go over"<<std::endl;
        return;
    }
    for( uint32_t i = 0 ; i < num ;  )
    {
        QString line = stream.readLine();
        QTextStream line_stream(&line);
        if(line==tr("=start="))
        {
            config.reset(new Config());
            config->add("Close_On_Finish","1");
        }else if(line==tr("=end=")){
            config->reload(configtxt);
            config_->updateBy(config);
            std::cerr<<"Starting Case "<<i+1<<"/"<<num<<std::endl;
            act->trigger();
            QApplication::processEvents();
            while( edit_thread_ || edit_widget_ )
            {
                QApplication::processEvents();
            }
            ui->mdiArea->closeAllSubWindows();
            mesh_views_.clear();
            inputs_.clear();
            labels_.clear();
            ++i;
        }else{
            QString str0;
            line_stream >> str0;
            if(str0.startsWith("#"))
            {
                continue;
            }
            QString str1 = line.remove(0,str0.size());
            str1 = str1.simplified();
            config_->add(str0.toStdString(),str1.toStdString());
        }
    }
}

QAction* MainWindow::getActionByText(const QString& txt)
{
    QList<QAction*> list;
    list = ui->menuEdit->actions();
    foreach(QAction* act,list)
    {
        if(act->text()==txt)return act;
    }
    list = ui->menuJRCS->actions();
    foreach(QAction* act,list)
    {
        if(act->text()==txt)return act;
    }
    list = ui->menuOptimization->actions();
    foreach(QAction* act,list)
    {
        if(act->text()==txt)return act;
    }
    return nullptr;
}

MainWindow::~MainWindow()
{
    delete ui;
}
