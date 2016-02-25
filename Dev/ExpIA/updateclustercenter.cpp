#include "updateclustercenter.h"
#include "extractpatchfeature.h"
#include "extractmesh.hpp"
#include "objectmodel.h"
#include <QTime>
#include "nanoflann.hpp"
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
    MeshKDTreeInterface<DefaultMesh> pts(pm);
    nanoflann::KDTreeSingleIndexAdaptor<
            nanoflann::L2_Simple_Adaptor<float,MeshKDTreeInterface<DefaultMesh>>,
            MeshKDTreeInterface<DefaultMesh>,
            3,arma::uword>
            kdtree(3,pts,nanoflann::KDTreeSingleIndexAdaptorParams(5));
    kdtree.buildIndex();
    double v = 0.0;
    float* p = (float*)om.points();

    arma::uword search_idx;
    float search_dist;

//    arma::fmat o_v_mat((float*)om.points(),3,om.n_vertices(),false,true);
//    arma::fmat p_v_mat((float*)pm.points(),3,pm.n_vertices(),false,true);

    arma::fmat o_n_mat((float*)om.vertex_normals(),3,om.n_vertices(),false,true);
    arma::fmat p_n_mat((float*)pm.vertex_normals(),3,pm.n_vertices(),false,true);

    arma::Mat<uint8_t> o_c_mat((uint8_t*)om.vertex_colors(),3,om.n_vertices(),false,true);
    arma::Mat<uint8_t> p_c_mat((uint8_t*)pm.vertex_colors(),3,pm.n_vertices(),false,true);

    for( size_t p_i = 0 ; p_i < om.n_vertices() ; ++ p_i )
    {
        kdtree.knnSearch(&p[3*p_i],1,&search_idx,&search_dist);
        double ncos = arma::dot(o_n_mat.col(p_i),p_n_mat.col(search_idx));
        arma::fvec oc = arma::conv_to<arma::fvec>::from(o_c_mat.col(p_i));
        arma::fvec pc = arma::conv_to<arma::fvec>::from(p_c_mat.col(search_idx));
        double cdist = arma::norm(pc - oc);
        v += ncos / ( 1.0 + cdist / 10.0 ) / ( 1.0 + 10.0*std::sqrt(search_dist) ) ;
    }
    return v;
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
    arma::mat X = arma::solve(Sw,Sb);
    arma::svd(U,s,V,X);
    int dim = 2;
    if(config_->has("Feature_dim"))dim = config_->getInt("Feature_dim");
    std::cerr<<"X:"<<X.n_rows<<","<<X.n_cols<<std::endl;
    std::cerr<<"V:"<<V.n_rows<<","<<V.n_cols<<std::endl;
    std::cerr<<"Base:"<<feature_base_.n_rows<<","<<feature_base_.n_cols<<std::endl;
    feature_base_.cols( 1 , dim ) = V.cols( 0 , dim - 1 );
}
//update the projection matrix by LDA
void UpdateClusterCenter::update()
{
    compute_Sw();
    compute_Sb();
    std::cerr<<"compute_proj"<<std::endl;
    compute_proj();
}

void UpdateClusterCenter::run()
{
    QTime timer;
    timer.restart();
    emit message(QString("evaluate patches"),0);
    evaluate_patches();
    emit message(QString("updating project matrix"),0);
    update();
    QString msg;
    msg = msg.sprintf("Cluster Center is Updated after %u ms",timer.elapsed());
    emit message(msg,0);
}

