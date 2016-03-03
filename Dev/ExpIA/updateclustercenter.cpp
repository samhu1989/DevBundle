#include "updateclustercenter.h"
#include "extractpatchfeature.h"
#include "extractmesh.hpp"
#include "objectmodel.h"
#include <QTime>
#include "nanoflann.hpp"
#include "OpenBLAS/lapacke.h"
#include "OpenBLAS/lapacke_utils.h"
bool UpdateClusterCenter::configure(Config::Ptr config)
{
    config_ = config;
    raw_feature_dim = 0;
    if(feature_base_.empty()){
        std::cerr<<"Err:Empty Feature Base"<<std::endl;
        return false;
    }
    if(objects_.empty())
    {
        std::cerr<<"Err:Empty Object List"<<std::endl;
        return false;
    }
    return true;
}

void UpdateClusterCenter::evaluate()
{
    invalid_objects_.clear();
    MeshBundle<DefaultMesh>::PtrList::iterator iter;
    arma::uword label_value = 1;
    arma::uword label_max = 0;
    LabelList::iterator liter;
    for(liter=labels_.begin();liter!=labels_.end();++liter)
    {
        arma::uword max = arma::max(*liter);
        if( max > label_max ) label_max = max;
    }
    std::vector<double> score;
    score.reserve(labels_.size());
    while(label_value <= label_max)
    {
        score.clear();
        patch_feature_.clear();
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
                if(patch_feature_.empty())patch_feature_ = feature;
                else patch_feature_.insert_cols( patch_feature_.n_cols - 1 , feature );
            }
            ++lindex;
        }
        patch_score_ = arma::rowvec(score);
        std::cerr<<"obj-"<<label_value<<std::endl;
        std::cerr<<"score("<<arma::min(patch_score_)<<","<<arma::max(patch_score_)<<","<<std::endl;
        ObjModel& model = *objects_[label_value-1];
        std::cerr<<"expected score:"<<model.GeoM_->mesh_.n_vertices()<<std::endl;

//        std::cerr<<"score num:"<<patch_score_.n_cols<<std::endl;
//        std::cerr<<"feature num:"<<patch_feature_.n_cols<<std::endl;
        if(validate_object())
        {
            select_samples(0.1*(model.GeoM_->mesh_.n_vertices()));
            compute_mi();
            compute_Hbi();
            compute_Hwi();
        }else{
            std::cerr<<"obj-"<<label_value<<"is invalid"<<std::endl;
            invalid_objects_.push_back(label_value);
        }
        ++ label_value;
    }
    if(!invalid_objects_.empty())remove_invalid();
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
        v = R*v;
        v.each_col() += t;
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
        v += ncos / ( 1.0 + 10.0*std::sqrt(search_dist) ) ;
    }
    return v;
}

bool UpdateClusterCenter::validate_object(double th)
{
    arma::uvec index = arma::find( patch_score_ > th );
    if(index.size()<patch_score_.size())return false;
    return true;
}

void UpdateClusterCenter::remove_invalid()
{
    ;
}

void UpdateClusterCenter::select_samples(double th)
{
    if(th<0)
    {
        double th = arma::mean(patch_score_);
        th = th > 0 ? th : 0 ;
    }
    arma::uvec selected = arma::find( patch_score_ >= th );
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
    mi.back() /= arma::accu(patch_score_);
}

void UpdateClusterCenter::compute_Hbi()
{
    Hbi.emplace_back(patch_feature_.n_rows,arma::fill::zeros);
    Hbi.back() = std::sqrt(double(Ni.back()))*mi.back();
}

void UpdateClusterCenter::compute_Hwi()
{
    Hwi.emplace_back(patch_feature_.n_rows,patch_feature_.n_cols,arma::fill::zeros);
    Hwi.back() = patch_feature_;
    Hwi.back().each_col() -= mi.back();
}

void UpdateClusterCenter::compute_Hw()
{
    std::vector<arma::mat>::iterator Hwiter;
    Hw = Hwi.front();
    for( Hwiter = (Hwi.begin()+1) ;Hwiter!=Hwi.end();++Hwiter)
    {
        Hw = arma::join_rows(Hw,*Hwiter);
    }
}

void UpdateClusterCenter::compute_Hb()
{
    Hb = arma::mat(Hbi.front().size(),Hbi.size(),arma::fill::zeros);
    std::vector<arma::vec>::iterator Hbiter;
    size_t index = 0;
    for(Hbiter=Hbi.begin();Hbiter!=Hbi.end();++Hbiter)
    {
        Hb.col(index) = *Hbiter;
        ++index;
    }
}

void UpdateClusterCenter::compute_Sw()
{
    Sw = Hw*Hw.t();
}

void UpdateClusterCenter::compute_Sb()
{
    Sb = Hb*Hb.t();
}

void UpdateClusterCenter::compute_base()
{
    arma::mat V;
    arma::vec s;
//    Sw.save("./debug/Sw.arma",arma::raw_ascii);
//    Sb.save("./debug/Sb.arma",arma::raw_ascii);
    arma::mat X = arma::pinv(Sw) * Sb;
    arma::eig_sym(s,V,X);
    int dim = 2;
    if(config_->has("Feature_dim"))dim = config_->getInt("Feature_dim");
    arma::uvec index = arma::sort_index(s,"descend");
    feature_base_.cols( 1 , dim ) = V.cols(index.head(dim));
}

void UpdateClusterCenter::compute_base_gsvd()
{
    ;
}

void UpdateClusterCenter::compute_center()
{
    std::vector<arma::vec>::iterator miter;
    size_t N = feature_base_.n_cols - 1;
    arma::mat proj = feature_base_.cols( 1, N );
    feature_centers_ = arma::mat(N,mi.size());
    size_t index = 0;
    for(miter=mi.begin();miter!=mi.end();++miter)
    {
        feature_centers_.col(index) = (((*miter).t())*proj).t();
        ++ index;
    }
}
//update the projection matrix by LDA
void UpdateClusterCenter::update()
{
    std::cerr<<"Hw"<<std::endl;
    compute_Hw();
    std::cerr<<"Hb"<<std::endl;
    compute_Hb();
    if( Hw.n_rows > arma::rank(Hw) )
    {
        std::cerr<<"Base"<<std::endl;
        compute_base_gsvd();
    }else{
        std::cerr<<"Sw"<<std::endl;
        compute_Sw();
        std::cerr<<"Sb"<<std::endl;
        compute_Sb();
        std::cerr<<"Base"<<std::endl;
        compute_base();
    }
    std::cerr<<"Center"<<std::endl;
    compute_center();
}

void UpdateClusterCenter::run()
{
    QTime timer;
    timer.restart();
    emit message(QString("evaluating"),0);
    evaluate();
    emit message(QString("updating"),0);
    update();
    QString msg;
    msg = msg.sprintf("Cluster Center is Updated after %u ms",timer.elapsed());
    emit message(msg,0);
}

