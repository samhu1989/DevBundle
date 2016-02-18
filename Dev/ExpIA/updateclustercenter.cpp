#include "updateclustercenter.h"
#include "extractpatchfeature.h"
#include "extractmesh.hpp"
#include "objectmodel.h"
#include <QTime>
bool UpdateClusterCenter::configure(Config::Ptr config)
{
    config_ = config;
    raw_feature_dim = 0;
}

void UpdateClusterCenter::evaluate_patches()
{
    patch_feature_.clear();
    MeshBundle<DefaultMesh>::PtrList::iterator iter;
    arma::uword label_value = 1;
    arma::uword label_max;
    std::vector<double> score;
    while(label_value <= label_max)
    {
        score.clear();
        arma::uword lindex = 0;
        for(iter=inputs_.begin();iter!=inputs_.end();++iter)
        {
            MeshBundle<DefaultMesh>::Ptr input = *iter;
            arma::uvec extracted_label = arma::find( label_value == labels_[lindex] );
            DefaultMesh extracted_mesh;
            if(!extracted_label.is_empty())
            {
                extractMesh<DefaultMesh,DefaultMesh>(input->mesh_,extracted_label,extracted_mesh);
                arma::vec feature;
                extract_patch_feature(extracted_mesh,feature,config_);
                feature -= feature_base_.col(0); //centralize the feature
                if(!raw_feature_dim)raw_feature_dim = feature.size();
                score.push_back( evaluate_patch( label_value - 1 , lindex , extracted_mesh ) );
                patch_feature_.insert_cols(patch_feature_.n_cols-1,feature);
            }
            ++lindex;
        }
        patch_score_ = arma::rowvec(score);
        select_samples();
        compute_mi();
        compute_Si();
        ++ label_value;
    }
}

double UpdateClusterCenter::evaluate_patch(uint64_t oidx,uint64_t fidx,DefaultMesh& patch)
{
    ObjModel& model = *objects_[oidx];
    ObjModel::T::Ptr& T_ptr = model.GeoT_[fidx];
    arma::fmat v((float*)patch.points(),3,patch.n_vertices(),false,true);
    arma::fmat vn((float*)patch.vertex_normals(),3,patch.n_vertices(),false,true);
    if(T_ptr&&0!=T_ptr.use_count())
    {
        arma::fmat R(T_ptr->R,3,3,false,true);
        arma::fvec t(T_ptr->t,3,false,true);
        v = R*v+t;
        vn = R*vn;
    }
    return match_patch(model.GeoM_->mesh_,patch);
}

double UpdateClusterCenter::match_patch(
        DefaultMesh& om,
        DefaultMesh& pm
        )
{
    ;
}

void UpdateClusterCenter::select_samples()
{
    double th = arma::mean(patch_score_);
    arma::uvec selected = arma::find( patch_score_ > th );
    Ni.push_back(selected.size());
    patch_score_ = patch_score_.cols(selected);
    patch_feature_ = patch_feature_.cols(selected);
}

void UpdateClusterCenter::compute_mi()
{
    mi.emplace_back(patch_feature_.n_rows,arma::fill::zeros);
    arma::mat tmp;
    tmp = patch_feature_;
    tmp.each_row() %= patch_score_;
    mi.back() = arma::sum(tmp,1);
    mi.back() /= arma::sum(patch_score_);
}

void UpdateClusterCenter::compute_Si()
{
    Si.emplace_back(patch_feature_.n_rows,patch_feature_.n_rows,arma::fill::zeros);
    arma::mat tmp;
    tmp = patch_feature_;
    tmp.each_col() -= mi.back();
    Si.back() = tmp * ( tmp.t() );
}

void UpdateClusterCenter::compute_Sw()
{
    Sw.reset();
    std::vector<arma::mat>::iterator Siter;
    for(Siter=Si.begin();Siter!=Si.end();++Siter)
    {
        if(Sw.is_empty())Sw = *Siter;
        else Sw += *Siter;
    }
}

void UpdateClusterCenter::compute_Sb()
{
    Sb = arma::mat(raw_feature_dim,raw_feature_dim,arma::fill::zeros);
    std::vector<arma::vec>::iterator miter;
    uint32_t index = 0;
    for(miter=mi.begin();miter!=mi.end();++miter)
    {
        Sb += Ni[index] * (*miter) * ( (*miter).t() );
        ++index;
    }
}

void UpdateClusterCenter::compute_proj()
{
    arma::mat U,V;
    arma::vec s;
    arma::svd(U,s,V,arma::solve(Sw,Sb));
    int dim = 2;
    if(config_->has("Feature_dim"))dim = config_->getInt("Feature_dim");
    feature_base_ = arma::mat(raw_feature_dim,1+dim);
    feature_base_.cols(1,dim) = V.cols( 0 , dim - 1 );
}

void UpdateClusterCenter::update()
{
    compute_Sw();
    compute_Sb();
    compute_proj();
}

void UpdateClusterCenter::run()
{
    QTime timer;
    timer.restart();
    emit message(QString("evaluate_patches"),0);
    evaluate_patches();
    update();
    QString msg;
    msg = msg.sprintf("Cluster Center is Updated after %u ms",timer.elapsed());
    emit message(msg,0);
}

