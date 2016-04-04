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
}

bool JRCSView::configure(Config::Ptr config)
{
    return init(config);
}

bool JRCSView::init(Config::Ptr config)
{
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
    MatPtrLst wv,wn;
    CMatPtrLst wc;
    LCMatPtrLst wl;
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

    }


    return true;
}

void JRCSView::move_worker_to_thread( JRCSThread* jrcs_worker )
{
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
    jrcs_thread_->start(QThread::HighPriority);
}

void JRCSView::finished()
{
    QThread* th = qobject_cast<QThread*>(sender());
    if(th&&th==jrcs_thread_)
    {
        jrcs_thread_->deleteLater();
        jrcs_thread_ = NULL;
    }
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
