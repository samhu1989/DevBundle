#include "extractbackground.h"
#include "ui_extractbackground.h"
#include "MeshPairViewerWidget.h"
#include "segmentationcore.h"
#include "extractplanethread.h"
#include "SAC/sac_plane.h"
#include "nanoflann.hpp"
ExtractBackground::ExtractBackground(
        std::vector<QWidget*>&inputs,
        QTimer& gl_timer,
        std::vector<arma::uvec>&labels,
        QWidget *parent
        ):
    QFrame(parent),
    inputs_(inputs),
    gl_timer_(gl_timer),
    labels_(labels),
    ui(new Ui::ExtractBackground)
{
    ui->setupUi(this);
    connect(ui->reset,SIGNAL(clicked(bool)),this,SLOT(reset_extract()));
    connect(ui->extract_n_plane,SIGNAL(clicked(bool)),this,SLOT(reset_extract_planes()));
    connect(ui->extract_n_plane,SIGNAL(clicked(bool)),this,SLOT(start_extract_planes()));
    connect(ui->extract_points,SIGNAL(clicked(bool)),this,SLOT(extract_points()));
    connect(ui->unextract_points,SIGNAL(clicked(bool)),this,SLOT(unextract_points()));
    connect(ui->extract_plane,SIGNAL(clicked(bool)),this,SLOT(extract_plane()));
}

bool ExtractBackground::configure(Config::Ptr)
{
    if(inputs_.empty())return false;
    if(labels_.size()!=inputs_.size())return false;
    return true;
}

ExtractBackground::~ExtractBackground()
{
    delete ui;
}

void ExtractBackground::start_extract_planes()
{
    MeshPairViewerWidget* w = (MeshPairViewerWidget*)(inputs_[current_frame_]);
    ExtractPlaneThread* th = new ExtractPlaneThread(w->first_ptr(),labels_[current_frame_]);
    th->setPlaneNumber(ui->n_planes->value());
    th->setThreshold(ui->dist->value());
    connect(th,SIGNAL(finished()),this,SLOT(finish_extract_planes()));
    th->start(QThread::NormalPriority);
}

void ExtractBackground::finish_extract_planes()
{
    ExtractPlaneThread* th = qobject_cast<ExtractPlaneThread*>(sender());
    if(!th)std::cerr<<"called by unknown object"<<std::endl;
    else{
       th->deleteLater();
    }
    ++current_frame_;
    if(current_frame_<inputs_.size())
    {
        start_extract_planes();
    }
}

void ExtractBackground::extract_points()
{
    std::vector<QWidget*>::iterator iter;
    uint64_t index = 0;
    for(iter=inputs_.begin();iter!=inputs_.end();++iter)
    {
        MeshPairViewerWidget* w = (MeshPairViewerWidget*)(*iter);
        extract_points(w->first_ptr(),w->first_selected(),labels_[index]);
        w->first_selected().clear();
        QApplication::processEvents();
        ++index;
    }
}
using namespace nanoflann;
void ExtractBackground::extract_points(
        MeshBundle<DefaultMesh>::Ptr ptr,
        const std::vector<arma::uword>& indices,
        arma::uvec&label
        )
{
    float* pts = (float*)ptr->mesh_.points();
    MeshKDTreeInterface<DefaultMesh> points(ptr->mesh_);
    KDTreeSingleIndexAdaptor<
            L2_Simple_Adaptor<float,MeshKDTreeInterface<DefaultMesh>>,
            MeshKDTreeInterface<DefaultMesh>,
            3,arma::uword>
            kdtree(3,points,KDTreeSingleIndexAdaptorParams(9));
    kdtree.buildIndex();
    arma::uvec result;
    std::vector<arma::uword>::const_iterator iter;
    for(iter=indices.cbegin();iter!=indices.cend();++iter)
    {
        arma::uvec found(ui->k_neighbors->value());
        arma::fvec dists(ui->k_neighbors->value());
        kdtree.knnSearch(
                &pts[3*(*iter)],
                ui->k_neighbors->value(),
                found.memptr(),
                dists.memptr()
                );
        result.insert_rows(0,1);
        result(0) = *iter;
        result = arma::join_cols(result,found);
    }
    label(result).fill(0);
    ptr->custom_color_.fromlabel(label);
}

void ExtractBackground::unextract_points()
{
    std::vector<QWidget*>::iterator iter;
    uint64_t index = 0;
    for(iter=inputs_.begin();iter!=inputs_.end();++iter)
    {
        MeshPairViewerWidget* w = (MeshPairViewerWidget*)(*iter);
        unextract_points(w->first_ptr(),w->first_selected(),labels_[index]);
        w->first_selected().clear();
        QApplication::processEvents();
        ++index;
    }
}

void ExtractBackground::unextract_points(
        MeshBundle<DefaultMesh>::Ptr ptr,
        const std::vector<arma::uword>& indices,
        arma::uvec&label
        )
{
    float* pts = (float*)ptr->mesh_.points();
    MeshKDTreeInterface<DefaultMesh> points(ptr->mesh_);
    KDTreeSingleIndexAdaptor<
            L2_Simple_Adaptor<float,MeshKDTreeInterface<DefaultMesh>>,
            MeshKDTreeInterface<DefaultMesh>,
            3,arma::uword>
            kdtree(3,points,KDTreeSingleIndexAdaptorParams(9));
    kdtree.buildIndex();
    arma::uvec result;
    std::vector<arma::uword>::const_iterator iter;
    for(iter=indices.cbegin();iter!=indices.cend();++iter)
    {
        arma::uvec found(ui->k_neighbors->value());
        arma::fvec dists(ui->k_neighbors->value());
        kdtree.knnSearch(
                &pts[3*(*iter)],
                ui->k_neighbors->value(),
                found.memptr(),
                dists.memptr()
                );
        result.insert_rows(0,1);
        result(0) = *iter;
        result = arma::join_cols(result,found);
    }
    label(result).fill(1);
    ptr->custom_color_.fromlabel(label);
}

void ExtractBackground::extract_plane()
{
    std::vector<QWidget*>::iterator iter;
    uint64_t index = 0;
    for(iter=inputs_.begin();iter!=inputs_.end();++iter)
    {
        MeshPairViewerWidget* w = (MeshPairViewerWidget*)(*iter);
        if(w->first_selected().size()!=3){
            ++index;
            continue;
        }
        extract_plane(w->first_ptr(),w->first_selected(),labels_[index]);
        w->first_selected().clear();
        QApplication::processEvents();
        ++index;
    }
}

void ExtractBackground::extract_plane(
        MeshBundle<DefaultMesh>::Ptr ptr,
        const std::vector<arma::uword>& indices,
        arma::uvec&label
        )
{
    Segmentation::SAC_Plane plane_model;
    arma::fmat input((float*)ptr->mesh_.points(),3,ptr->mesh_.n_vertices(),false,true);
    plane_model.input(input);
    arma::uvec inliers(indices);
    arma::vec coeff;
    plane_model.computeModel(inliers,coeff);
    arma::uvec result;
    plane_model.selectWithinDistance(coeff,ui->dist->value(),result);
    label(result).fill(0);
    ptr->custom_color_.fromlabel(label);
}

void ExtractBackground::reset_extract()
{
    std::vector<QWidget*>::iterator iter;
    uint64_t index = 0;
    for(iter=inputs_.begin();iter!=inputs_.end();++iter)
    {
        MeshPairViewerWidget* w = (MeshPairViewerWidget*)(*iter);
        labels_[index].fill(0);
        w->first_ptr()->custom_color_.fromlabel(labels_[index]);
        ++index;
    }
}
