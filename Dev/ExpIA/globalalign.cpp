#include "globalalign.h"
#include "ui_globalalign.h"
#include "filter.h"
#include "registrationcore.h"
GlobalAlign::GlobalAlign(GlobalAlign::IMeshList& input,QWidget* parent) :
    QFrame(parent),
    inputs_(input),
    geo_thread_(NULL),
    ui(new Ui::GlobalAlign)
{
    ui->setupUi(this);
    geo_view_ = new MeshListViewerWidget();
    geo_view_->setMinimumSize(640,480);
    ui->verticalLayout->addWidget(geo_view_);

    connect(ui->load_down_sampled,SIGNAL(clicked(bool)),this,SLOT(loadDownSampled()));
    connect(ui->align_each_other,SIGNAL(clicked(bool)),this,SLOT(start_AlignEachOther()));
    connect(ui->align_up,SIGNAL(clicked(bool)),this,SLOT(alignUptoZ()));

    gl_timer.setSingleShot(false);
    connect(&gl_timer,SIGNAL(timeout()),geo_view_,SLOT(updateGL()));
    gl_timer.start(100);
}

bool GlobalAlign::configure(Config::Ptr config)
{
    config_ = config;
    if(inputs_.empty())return false;
    return true;
}

void GlobalAlign::loadDownSampled(void)
{
    IMeshList::iterator iter;
    for(iter=inputs_.begin();iter!=inputs_.end();++iter)
    {
        MeshBundle<DefaultMesh>::Ptr ptr;
        ptr = *iter;
        MeshBundle<DefaultMesh>::Ptr downSampled(new MeshBundle<DefaultMesh>);
        downSampled->mesh_.request_vertex_colors();
        downSampled->mesh_ = ptr->mesh_;
        Filter::OctreeGrid<DefaultMesh> filter;
        filter.set_seed_resolution(ui->resolution_spin_box->value());
        filter.extract(downSampled->mesh_);
        geo_view_->list().push_back(downSampled);
    }
    geo_view_->set_center_at_mesh(geo_view_->list().back()->mesh_);
}

void GlobalAlign::start_AlignEachOther(void)
{
    if(geo_thread_)
    {
        QString msg = "Please Wait Till the End of Last Registration\n";
        QMessageBox::critical(this, windowTitle(), msg);
        return;
    }
    JRMPC_Thread* th = new JRMPC_Thread();
    if(!th->init(geo_view_->list(),config_))
    {
        QString msg = "Fail to Initialize the Registration:\n '";
        msg += QString::fromStdString(th->errorString());
        QMessageBox::critical( NULL, windowTitle(), msg);
        return;
    }else{
        geo_view_->reset_center();
        geo_view_->show_back();
        connect(th,SIGNAL(finished()),this,SLOT(finish_AlignEachOther()));
    }
    geo_thread_ = th;
    QString name;
    name = "Aligning";
    geo_thread_->setObjectName(name);
    emit message(name,0);
    geo_thread_->start(QThread::HighestPriority);
}

void GlobalAlign::finish_AlignEachOther(void)
{
    QThread* th = qobject_cast<QThread*>(sender());
    if(th)
    {
        QString msg = th->objectName() + " is Finished";
        if(th==geo_thread_)
        {
            while(geo_thread_->isRunning())
            {
                geo_thread_->terminate();
                QApplication::processEvents();
            }
            JRMPC_Thread* jrmpcth = qobject_cast<JRMPC_Thread*>(geo_thread_);
            Registration::JRMPC<DefaultMesh>::ResPtr r = jrmpcth->result();
            arma::uword index = 0;
            IMeshList::iterator iter;
            for(iter=inputs_.begin();iter!=inputs_.end();++iter)
            {
                arma::fmat R = *(r->Rs[index]);
                arma::fvec t = *(r->ts[index]);
                float* pdata = (float*)(*iter)->mesh_.points();
                arma::fmat pts(pdata,3,(*iter)->mesh_.n_vertices(),false,true);
                pts = R*pts;
                pts.each_col() += t;
                ++index;
            }
            geo_thread_->deleteLater();
            geo_thread_ = NULL;
        }
        emit message(msg,0);
    }
}

void GlobalAlign::alignUptoZ(void)
{
    if(geo_view_->current_visible_num()!=1)return;
    if( 4 != geo_view_->current_selected().size() )return;
    arma::uvec selected(geo_view_->current_selected());
    float* pdata = (float*)geo_view_->list()[geo_view_->current_mesh_start()]->mesh_.points();
    uint32_t N = geo_view_->list()[geo_view_->current_mesh_start()]->mesh_.n_vertices();
    arma::fmat points(pdata,3,N,false,true);
    arma::fmat picked = points.cols(selected);
    arma::fvec t = arma::mean(picked.cols(0,2),1);
    arma::fvec zUp(3,arma::fill::zeros);
    zUp(2) = 1.0;
    arma::fvec norm = arma::cross( picked.col(1) - picked.col(0) , picked.col(2) - picked.col(0) );
    arma::fvec tmp = picked.col(3) - t;
    if( std::abs(arma::dot( norm , tmp )) < 0 )norm*=-1.0;
    norm = arma::normalise(norm);
    arma::fmat R(3,3);
    getRotationFromTwoUnitVectors(zUp,norm,R);
    //transform downsampled
    IMeshList::iterator iter;
    IMeshList::iterator b = geo_view_->list().begin();
    IMeshList::iterator e = geo_view_->list().end();
    e--;//avoid last one
    for(iter=b;iter!=e;++iter)
    {
        float* pdata = (float*)(*iter)->mesh_.points();
        arma::fmat pts(pdata,3,(*iter)->mesh_.n_vertices(),false,true);
        pts.each_col() -= t;
        pts = R*pts;
    }
    //transform inputs
    for(iter=inputs_.begin();iter!=inputs_.end();++iter)
    {
        float* pdata = (float*)(*iter)->mesh_.points();
        arma::fmat pts(pdata,3,(*iter)->mesh_.n_vertices(),false,true);
        pts.each_col() -= t;
        pts = R*pts;
    }
    geo_view_->current_selected().clear();
}

GlobalAlign::~GlobalAlign()
{
    delete ui;
}
