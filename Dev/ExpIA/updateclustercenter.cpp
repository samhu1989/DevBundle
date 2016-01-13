#include "updateclustercenter.h"
#include "extractpatchfeature.h"
#include "extractmesh.hpp"
#include "objectmodel.h"
#include <QTime>
bool UpdateClusterCenter::configure(Config::Ptr config)
{
    config_ = config;
}

void UpdateClusterCenter::evaluate_patches()
{
    patch_frame_.clear();
    patch_class_.clear();
    patch_feature_.clear();
    MeshBundle<DefaultMesh>::PtrList::iterator iter;
    size_t lindex = 0;
    int feature_dim;
    for(iter=inputs_.begin();iter!=inputs_.end();++iter)
    {
        MeshBundle<DefaultMesh>::Ptr input = *iter;
        arma::uword label_value = 1;
        arma::uword label_max = arma::max(labels_[lindex]);
        while( label_value <= label_max )
        {
            arma::uvec extracted_label = arma::find( label_value == labels_[lindex] );
            DefaultMesh extracted_mesh;
            if(!extracted_label.is_empty())
            {
                extractMesh<DefaultMesh,DefaultMesh>(input->mesh_,extracted_label,extracted_mesh);
                arma::vec feature;
                extract_patch_feature(extracted_mesh,feature,config_);
                if(!raw_feature_dim)raw_feature_dim = feature.size();
                patch_frame_.push_back( lindex );
                patch_class_.push_back( label_value );
                patch_score_.push_back( evaluate_patch( label_value - 1 , lindex , extracted_mesh ) );
                patch_feature_.insert_cols(patch_feature_.n_cols-1,patch_feature_);
            }
            ++ label_value;
        }
        ++lindex;
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

void UpdateClusterCenter::evaluate_objects()
{
    ;
}

void UpdateClusterCenter::update_proj()
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
    update_proj();
}

void UpdateClusterCenter::run()
{
    QTime timer;
    timer.restart();
    emit message(QString("evaluate_patches"),0);
    evaluate_patches();
    emit message(QString("evaluate_objects"),0);
    evaluate_objects();
    update();
    QString msg;
    msg = msg.sprintf("Cluster Center is Updated after %u ms",timer.elapsed());
    emit message(msg,0);
}

