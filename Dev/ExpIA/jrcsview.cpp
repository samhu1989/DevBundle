#include "jrcsview.h"
#include "ui_jrcsview.h"
#include <QThread>
#include <armadillo>
JRCSView::JRCSView(
        MeshList& inputs,
        LabelList& labels,
        ModelList& objects,
        QWidget *parent
        ):
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
        return false;
    }
    input(jrcs_worker_);
    if(!allocate_x(jrcs_worker_))return false;
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

        std::cerr<<"face num:"<<mesh.n_faces()<<std::endl;

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

    std::cerr<<"face num:"<<mesh.n_faces()<<std::endl;

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
        jrcs_worker_->deleteLater();
        jrcs_worker_ = NULL;
        jrcs_thread_->deleteLater();
        jrcs_thread_ = NULL;
    }
    std::cout<<msg.toStdString()<<std::endl;
    emit message(msg,1000);
    connect(ui->spinBox,SIGNAL(valueChanged(int)),this,SLOT(align(int)));
    ui->spinBox->show();
    move(pos().x()+1,pos().y());
    align(ui->spinBox->value());
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

JRCSView::~JRCSView()
{
    geo_view_->close();
    ui->gridLayout->removeWidget(geo_view_);
    geo_view_->deleteLater();
    delete ui;
}
