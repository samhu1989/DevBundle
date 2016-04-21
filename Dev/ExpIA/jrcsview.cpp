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
    ui->gridLayout->addWidget(geo_view_);
    jrcs_thread_ = NULL;
    connect(&t_,SIGNAL(timeout()),geo_view_,SLOT(updateGL()));
    t_.setSingleShot(false);
}

bool JRCSView::configure(Config::Ptr config)
{
    return init(config);
}

bool JRCSView::init(Config::Ptr config)
{
    std::cerr<<"JRCSView::init"<<std::endl;
    JRCSThread* worker = new JRCSThread();
    if(!worker->configure(config))
    {
        return false;
    }
    input(worker);
    if(!allocate_x(worker))return false;
    move_worker_to_thread(worker);
    return true;
}

void JRCSView::input( JRCSThread* jrcs_worker_ )
{
    std::cerr<<"JRCSView::input"<<std::endl;
    MeshList::iterator iter;
    MatPtrLst vv,vn;
    CMatPtrLst vc;
    LCMatPtrLst vl;
    for(iter=inputs_.begin();iter!=inputs_.end();++iter)
    {
        MeshBundle<DefaultMesh>& mesh = **iter;
        int N = mesh.mesh_.n_vertices();
        vv.emplace_back(new arma::fmat((float*)mesh.mesh_.points(),3,N,false,true));
        vn.emplace_back(new arma::fmat((float*)mesh.mesh_.vertex_normals(),3,N,false,true));
        vc.emplace_back(new arma::Mat<uint8_t>((uint8_t*)mesh.mesh_.vertex_colors(),3,N,false,true));
        vl.emplace_back(new arma::Col<uint32_t>((uint32_t*)mesh.custom_color_.vertex_colors(),N,false,true));
    }
    jrcs_worker_->input(vv,vn,vc,vl);
}

bool JRCSView::allocate_x( JRCSThread* jrcs_worker_ )
{
    std::cerr<<"JRCSView::allocate_x"<<std::endl;
    MatPtrLst wv,wn;
    CMatPtrLst wc;
    int k = jrcs_worker_->get_k();
    if(k<0)return false;
    int M = inputs_.size();
    for(int i = 0 ; i < M ; i++ )
    {
        geo_view_->list().emplace_back(new MeshBundle<DefaultMesh>());
        DefaultMesh& mesh = geo_view_->list().back()->mesh_;
        for( int j=0 ; j < k ; ++j )
        {
            mesh.add_vertex(DefaultMesh::Point(0,0,0));
        }
        mesh.request_vertex_normals();
        mesh.request_vertex_colors();
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
    std::cerr<<"k:"<<k<<std::endl;
    for( int j=0 ; j < k ; ++j )
    {
        mesh.add_vertex(DefaultMesh::Point(0,0,0));
    }
    mesh.request_vertex_normals();
    mesh.request_vertex_colors();
    xv = std::make_shared<arma::fmat>((float*)mesh.points(),3,mesh.n_vertices(),false,true);
    xn = std::make_shared<arma::fmat>((float*)mesh.vertex_normals(),3,mesh.n_vertices(),false,true);
    xc = std::make_shared<arma::Mat<uint8_t>>((uint8_t*)mesh.vertex_colors(),3,mesh.n_vertices(),false,true);

    std::cerr<<"JRCSView::allocate_x:reseting x"<<std::endl;
    jrcs_worker_->resetx(xv,xn,xc);
    geo_view_->show_back();

    return true;
}

void JRCSView::move_worker_to_thread( JRCSThread* jrcs_worker )
{
    std::cerr<<"JRCSView::move_worker_to_thread"<<std::endl;
    jrcs_thread_ = new QThread();
    jrcs_worker->moveToThread(jrcs_thread_);
    connect(jrcs_thread_,SIGNAL(started()),jrcs_worker,SLOT(process()));
    connect(jrcs_worker,SIGNAL(end()),jrcs_worker,SLOT(deleteLater()));
    connect(jrcs_worker,SIGNAL(destroyed(QObject*)),jrcs_thread_,SLOT(quit()));
    connect(jrcs_worker,SIGNAL(message(QString,int)),this,SLOT(passMessage(QString,int)));
    connect(jrcs_thread_,SIGNAL(finished()),this,SLOT(finished()));
}

void JRCSView::start()
{
    std::cerr<<"JRCSView::start"<<std::endl;
    jrcs_thread_->start(QThread::HighPriority);
    t_.start(300);
}

void JRCSView::finished()
{
    std::cerr<<"JRCSView::finished"<<std::endl;
    QThread* th = qobject_cast<QThread*>(sender());
    if(th&&th==jrcs_thread_)
    {
        jrcs_thread_->deleteLater();
        jrcs_thread_ = NULL;
    }
    emit message(tr("JRCSView is finished"),0);
    emit closeInMdi(parentWidget());
}

void JRCSView::passMessage(QString msg,int t)
{
    emit message(msg,t);
}

JRCSView::~JRCSView()
{
    geo_view_->close();
    ui->gridLayout->removeWidget(geo_view_);
    geo_view_->deleteLater();
    delete ui;
}
