#include "jrcsview.h"
#include "ui_jrcsview.h"
#include <QThread>
#include <armadillo>
#include <iomanip>
JRCSView::JRCSView(
        MeshList& inputs,
        LabelList& labels,
        ModelList& objects,
        QWidget *parent
        ):
    close_on_finish_(false),
    inputs_(inputs),
    labels_(labels),
    objects_(objects),
    QFrame(parent),
    ui(new Ui::JRCSView)
{
    ui->setupUi(this);
    geo_view_ = new MeshListViewerWidget();
    geo_view_->setMinimumSize(640,480);
    ui->gridLayout_2->addWidget(geo_view_);
    jrcs_thread_ = NULL;
    jrcs_worker_ = new JRCSThread();
    connect(&t_,SIGNAL(timeout()),geo_view_,SLOT(updateGL()));
    t_.setSingleShot(false);
}

bool JRCSView::configure(Config::Ptr config)
{
    std::cerr<<"JRCSView::configure"<<std::endl;
    if(config->has("Close_On_Finish"))
    {
        if(config->getInt("Close_On_Finish"))
        {
            close_on_finish_ = true;
        }
    }
    return init(config);
}

bool JRCSView::init(Config::Ptr config)
{
    std::cerr<<"JRCSView::init"<<std::endl;
    std::cerr<<"On "<<jrcs_worker_->get_method_name()<<std::endl;
    ui->spinBox->hide();
    connect(jrcs_worker_,SIGNAL(message(QString,int)),this,SLOT(passMessage(QString,int)));
    connect(&t_,SIGNAL(timeout()),jrcs_worker_,SLOT(get_iter_info()));
    if(!jrcs_worker_->configure(config))
    {
        std::cerr<<"Failed to Init "<<jrcs_worker_->get_method_name()<<std::endl;
        return false;
    }
    input(jrcs_worker_);
    if(!allocate_x(jrcs_worker_)){
        std::cerr<<"Failed to allocate_x for "<<jrcs_worker_->get_method_name()<<std::endl;
        return false;
    }
    std::cerr<<"---"<<std::endl;
    move_worker_to_thread(jrcs_worker_);
    return true;
}


void JRCSView::input( JRCSThread* jrcs_worker_ )
{
    std::cerr<<"JRCSView::input"<<std::endl;
    MeshList::iterator iter;
    MatPtrLst vv,vn;
    CMatPtrLst vc;
    LCMatPtrLst vlc;
    bool use_init_label_ = true;
    uint32_t index = 0;
    for(iter=inputs_.begin();iter!=inputs_.end();++iter)
    {
        MeshBundle<DefaultMesh>& mesh = **iter;
        int N = mesh.mesh_.n_vertices();
        vv.emplace_back(new arma::fmat((float*)mesh.mesh_.points(),3,N,false,true));
        vn.emplace_back(new arma::fmat((float*)mesh.mesh_.vertex_normals(),3,N,false,true));
        vc.emplace_back(new arma::Mat<uint8_t>((uint8_t*)mesh.mesh_.vertex_colors(),3,N,false,true));
        vlc.emplace_back(new arma::Col<uint32_t>((uint32_t*)mesh.custom_color_.vertex_colors(),N,false,true));
        arma::uvec& lbl = labels_[index];
        if( arma::max(lbl) == arma::min(lbl) ) use_init_label_ = false;
        ++ index;
    }
    if(use_init_label_)
    {
        std::cerr<<"using init labels"<<std::endl;
        LMatPtrLst vl;
        std::vector<arma::uvec>::iterator liter;
        for(liter=labels_.begin();liter!=labels_.end();++liter)
        {
            vl.emplace_back(new arma::uvec(*liter));
        }
        jrcs_worker_->input_with_label(vv,vn,vc,vlc,vl);
    }else jrcs_worker_->input(vv,vn,vc,vlc);
}


void build_mesh(DefaultMesh& mesh,const int k,const int obj_n)
{
    std::vector<DefaultMesh::VertexHandle>  face_vhandles_a,face_vhandles_b;
    face_vhandles_a.reserve(3);
    face_vhandles_a.reserve(3);
    DefaultMesh::Point p(std::numeric_limits<float>::quiet_NaN(),std::numeric_limits<float>::quiet_NaN(),std::numeric_limits<float>::quiet_NaN());
    for( int j=0 ; j < k ;  )
    {
        if( (j + 4) <= 5*4*obj_n )
        {
            face_vhandles_a.push_back(mesh.add_vertex(p));
            face_vhandles_a.push_back(mesh.add_vertex(p));
            face_vhandles_a.push_back(mesh.add_vertex(p));
            face_vhandles_b.push_back(mesh.add_vertex(p));
            face_vhandles_b.push_back(face_vhandles_a[0]);
            face_vhandles_b.push_back(face_vhandles_a[2]);
            mesh.add_face(face_vhandles_a);
            mesh.add_face(face_vhandles_b);
            face_vhandles_a.clear();
            face_vhandles_b.clear();
            j += 4;
        }else{
            mesh.add_vertex(p);
            ++ j;
        }
    }

    mesh.request_face_normals();
    mesh.request_face_colors();
    mesh.request_vertex_normals();
    mesh.request_vertex_colors();
}

bool JRCSView::allocate_x( JRCSThread* jrcs_worker_ )
{
    std::cerr<<"JRCSView::allocate_x"<<std::endl;
    MatPtrLst wv,wn;
    CMatPtrLst wc;
    int k = jrcs_worker_->get_k();
    if(k<0){
        std::cerr<<"JRCSView::invalid k:"<<k<<std::endl;
        return false;
    }
    int M = inputs_.size();
    for(int i = 0 ; i < M ; i++ )
    {
        geo_view_->list().emplace_back(new MeshBundle<DefaultMesh>());
        DefaultMesh& mesh = geo_view_->list().back()->mesh_;

        build_mesh(mesh,k,jrcs_worker_->get_obj_n());

//        std::cerr<<"face num:"<<mesh.n_faces()<<std::endl;

        wv.emplace_back(new arma::fmat((float*)mesh.points(),3,mesh.n_vertices(),false,true));
        wn.emplace_back(new arma::fmat((float*)mesh.vertex_normals(),3,mesh.n_vertices(),false,true));
        wc.emplace_back(new arma::Mat<uint8_t>((uint8_t*)mesh.vertex_colors(),3,mesh.n_vertices(),false,true));
    }
    jrcs_worker_->resetw(wv,wn,wc);
    std::cerr<<"JRCSView::allocate_x:donewv"<<std::endl;
    MatPtr xv,xn;
    CMatPtr xc;
    geo_view_->list().emplace_back(new MeshBundle<DefaultMesh>());
    DefaultMesh& mesh = geo_view_->list().back()->mesh_;

    build_mesh(mesh,k,jrcs_worker_->get_obj_n());

//    std::cerr<<"face num:"<<mesh.n_faces()<<std::endl;

    xv = std::make_shared<arma::fmat>((float*)mesh.points(),3,mesh.n_vertices(),false,true);
    xn = std::make_shared<arma::fmat>((float*)mesh.vertex_normals(),3,mesh.n_vertices(),false,true);
    xc = std::make_shared<arma::Mat<uint8_t>>((uint8_t*)mesh.vertex_colors(),3,mesh.n_vertices(),false,true);

    std::cerr<<"JRCSView::allocate_x:reseting x"<<std::endl;
    jrcs_worker_->resetx(xv,xn,xc);
    mesh.update_face_normals();
    std::cerr<<"JRCSView::allocate_x:done reseting x"<<std::endl;
    geo_view_->reset_center();
    geo_view_->show_back();

    return true;
}

void JRCSView::move_worker_to_thread( JRCSThread* jrcs_worker )
{
    emit message(tr("JRCSView::move_worker_to_thread"),1000);
    jrcs_thread_ = new QThread();
    jrcs_worker->moveToThread(jrcs_thread_);
    connect(jrcs_thread_,SIGNAL(started()),jrcs_worker,SLOT(process()));
    connect(jrcs_worker,SIGNAL(end()),jrcs_thread_,SLOT(quit()));
    connect(jrcs_thread_,SIGNAL(finished()),this,SLOT(finished()));
}

void JRCSView::start()
{
    emit message(tr("Starting"),1000);
    if(!jrcs_worker_->input_extra(inputs_))
    {
        emit message(tr("Failed at input extra: perhaps missing some inputs"),0);
        std::cerr<<"Failed at input extra: perhaps missing some inputs"<<std::endl;
    }
    jrcs_thread_->start(QThread::NormalPriority);
    t_.start(751);
    time.restart();
}

void JRCSView::finished()
{
    std::cerr<<"JRCSView::finished"<<std::endl;
    QThread* th = qobject_cast<QThread*>(sender());
    QString msg;
    msg = msg.sprintf("%s is finished: %d iterations in %u ms",jrcs_worker_->get_method_name().c_str(),jrcs_worker_->get_iter(),time.elapsed());
    if(th&&th==jrcs_thread_)
    {
        t_.disconnect(geo_view_,SLOT(updateGL()));
        t_.stop();
        ui->spinBox->setMinimum(0);
        ui->spinBox->setMaximum(jrcs_worker_->get_obj_n()-1);
        ui->spinBox->setReadOnly(false);
        jrcs_worker_->get_lbl(labels_);
        jrcs_worker_->get_rt(rt_);
        jrcs_worker_->get_order(order_);
        jrcs_worker_->deleteLater();
        jrcs_worker_ = NULL;
        jrcs_thread_->deleteLater();
        jrcs_thread_ = NULL;
    }
    std::cout<<msg.toStdString()<<std::endl;
    emit message(msg,1000);
    if(close_on_finish_)
    {
        emit closeInMdi(this);
    }

    save_rt();
    if(!order_.empty())save_order();
    save_centroids();
    move(pos().x()+1,pos().y());
    ui->spinBox->show();
    connect(ui->spinBox,SIGNAL(valueChanged(int)),this,SLOT(align(int)));
    align(ui->spinBox->value());
    this->updateGeometry();
    this->update();
}

void JRCSView::passMessage(QString msg,int t)
{
    emit message(msg,t);
}

void JRCSView::align(int i)
{
    geo_view_->list().clear();
    MeshList::iterator iter;
    int idx = 0;
    for(iter=inputs_.begin();iter!=inputs_.end();++iter)
    {
        geo_view_->list().emplace_back(new MeshBundle<DefaultMesh>());
        geo_view_->list().back()->mesh_;
        geo_view_->list().back()->mesh_.request_vertex_colors();
        geo_view_->list().back()->mesh_.request_vertex_normals();
        geo_view_->list().back()->mesh_ = (*iter)->mesh_;
        DefaultMesh& mesh = geo_view_->list().back()->mesh_;
        arma::fmat R(rt_[idx][i].R,3,3,false,true);
        arma::fvec t(rt_[idx][i].t,3,false,true);
        arma::fmat vv((float*)mesh.points(),3,mesh.n_vertices(),false,true);
        arma::fmat vn((float*)mesh.vertex_normals(),3,mesh.n_vertices(),false,true);
        vv.each_col() -= t;
        vv = R.i()*vv;
        vn = R.i()*vn;
        ++idx;
    }
    geo_view_->show_back();
}

void JRCSView::save_rt()
{
    QString dirName;
    while(dirName.isEmpty())
    {
        dirName = QFileDialog::getExistingDirectory(
                this,
                tr("Save RT"),
                tr("../Dev_Data/")
                );
    }
    QDir dir;
    dir.setPath(dirName);
    for(int idx=0;idx<rt_.size();++idx)
    {
        QString fileName = dir.absoluteFilePath(
                    QString::fromStdString(inputs_[idx]->name_+".txt")
                    );
        QString fileName2 = dir.absoluteFilePath(
                    QString::fromStdString(inputs_[idx]->name_+".FromX.txt")
                    );
        std::ofstream out;
        std::ofstream out2;
        out.open(fileName.toStdString());
        out2.open(fileName2.toStdString());
        for(int o=0;o<rt_[idx].size();++o)
        {
            arma::fmat Rout(3,3,arma::fill::eye);
            arma::fvec tout(3,arma::fill::zeros);
            if(idx>0)
            {
                Rout = arma::fmat(rt_[idx-1][o].R,3,3,true,true);
                tout = arma::fmat(rt_[idx-1][o].t,3,true,true);
            }else{
                Rout = arma::fmat(rt_.back()[o].R,3,3,true,true);
                tout = arma::fmat(rt_.back()[o].t,3,true,true);
            }
            arma::fmat R(rt_[idx][o].R,3,3,false,true);
            arma::fvec t(rt_[idx][o].t,3,false,true);

            Rout = R * Rout.i();
            tout = Rout*(- tout)+t;
            out << Rout(0,0)<<" "<< Rout(0,1)<<" "<<Rout(0,2)<<" "<<tout(0)<<std::endl;
            out << Rout(1,0)<<" "<< Rout(1,1)<<" "<<Rout(1,2)<<" "<<tout(1)<<std::endl;
            out << Rout(2,0)<<" "<< Rout(2,1)<<" "<<Rout(2,2)<<" "<<tout(2)<<std::endl;
            out2<<std::setprecision(8)<<std::scientific;
            out2 << R(0,0)<<" "<< R(0,1)<<" "<<R(0,2)<<" "<<t(0)<<std::endl;
            out2 << R(1,0)<<" "<< R(1,1)<<" "<<R(1,2)<<" "<<t(1)<<std::endl;
            out2 << R(2,0)<<" "<< R(2,1)<<" "<<R(2,2)<<" "<<t(2)<<std::endl;
        }
        out.close();
        out2.close();
    }

}

void JRCSView::save_centroids()
{
    QString fileName;
    while(fileName.isEmpty())
    {
        fileName = QFileDialog::getSaveFileName(
                this,
                tr("Save Centroids"),
                tr("../Dev_Data/"),
                tr(
                    "PLY Files (*.ply);;"
                    "OBJ Files (*.obj);;"
                    "OFF Files (*.off);;"
                    "STL Files (*.stl);;"
                    "All Files (*)"
                   )
                );
    }
    DefaultMesh& m = geo_view_->list().back()->mesh_;
    OpenMesh::IO::Options opt;
    opt+=OpenMesh::IO::Options::Binary;
    if(m.has_vertex_colors())opt+=OpenMesh::IO::Options::VertexColor;
    if(m.has_vertex_normals())opt+=OpenMesh::IO::Options::VertexNormal;
    if(!OpenMesh::IO::write_mesh(m,fileName.toStdString(),opt,13)){
        std::cerr<<"can't save to:"<<fileName.toStdString()<<std::endl;
        return;
    }

}

void JRCSView::save_order()
{
    QString dirName;
    while(dirName.isEmpty())
    {
        dirName = QFileDialog::getExistingDirectory(
                this,
                tr("Save Order"),
                tr("../Dev_Data/")
                );
    }
    std::vector<arma::uvec>::iterator iter;
    std::vector<MeshBundle<DefaultMesh>::Ptr>::iterator miter;
    QDir dir;
    dir.setPath(dirName);
    miter = inputs_.begin();
    iter = order_.begin();
    QString filepath = dir.absoluteFilePath(
                QString::fromStdString("X.order.arma")
                );
    if(!iter->save(filepath.toStdString(),arma::arma_binary))
    {
        QString msg = "Failed to Save "+filepath+"\n";
        QMessageBox::critical(this, windowTitle(), msg);
        return;
    }
    ++iter;
    for( ;iter!=order_.end();++iter)
    {
        QString filepath = dir.absoluteFilePath(
                    QString::fromStdString((*miter)->name_+".order.arma")
                    );
        if(!iter->save(filepath.toStdString(),arma::arma_binary))
        {
            QString msg = "Failed to Save "+filepath+"\n";
            QMessageBox::critical(this, windowTitle(), msg);
            return;
        }
        ++miter;
        if(miter==inputs_.end())break;
    }
}

JRCSView::~JRCSView()
{
    if(jrcs_thread_)
    {
        if(jrcs_thread_->isRunning())
        {
            jrcs_thread_->terminate();
            jrcs_thread_->exit(-1);
            jrcs_thread_->quit();
        }
        jrcs_thread_->deleteLater();
    }
    if(jrcs_worker_)
    {
        delete jrcs_worker_;
    }
    geo_view_->close();
    ui->gridLayout->removeWidget(geo_view_);
    geo_view_->deleteLater();
    delete ui;
}
