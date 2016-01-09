#include "extractpatchfeature.h"
#include "nanoflann.hpp"
void extract_patch_feature(DefaultMesh&mesh,arma::vec&feature,Config::Ptr config_)
{
    arma::vec block_feature;//length width height
    Feature::BlockBasedFeature<DefaultMesh> block_feature_extractor;
    block_feature_extractor.extract(mesh,block_feature);
    arma::vec color_hist;
    if("RGB"==config_->getString("Color_Space"))
    {
        Feature::ColorHistogramRGB<DefaultMesh>
                color_hist_extractor(
                    config_->getInt("RGB_r_bin_num"),
                    config_->getInt("RGB_g_bin_num"),
                    config_->getInt("RGB_b_bin_num")
                    );
        color_hist_extractor.extract(mesh,color_hist);
    }else if("Lab"==config_->getString("Color_Space")){
        Feature::ColorHistogramLab<DefaultMesh>
                color_hist_extractor(
                    config_->getInt("Lab_L_bin_num"),
                    config_->getInt("Lab_a_bin_num"),
                    config_->getInt("Lab_b_bin_num")
                    );
        color_hist_extractor.extract(mesh,color_hist);
    }
    if(config_->has("Feature_height_w"))
    {
        block_feature(2) *= config_->getFloat("Feature_height_w");
//        block_feature(5) *= config_->getFloat("Feature_height_w");
//        block_feature(8) *= config_->getFloat("Feature_height_w");
//        block_feature(11) *= config_->getFloat("Feature_height_w");
    }
    if(config_->has("Feature_length_width_w"))
    {
        block_feature(0) *= config_->getFloat("Feature_length_width_w");
        block_feature(1) *= config_->getFloat("Feature_length_width_w");
//        block_feature(3) *= config_->getFloat("Feature_length_width_w");
//        block_feature(4) *= config_->getFloat("Feature_length_width_w");
//        block_feature(6) *= config_->getFloat("Feature_length_width_w");
//        block_feature(7) *= config_->getFloat("Feature_length_width_w");
//        block_feature(9) *= config_->getFloat("Feature_length_width_w");
//        block_feature(10) *= config_->getFloat("Feature_length_width_w");
    }
    if(config_->has("Feature_color_hist_w"))
    {
        color_hist *= config_->getFloat("Feature_color_hist_w");
    }
    feature = arma::join_cols(block_feature,color_hist);
//    feature = color_hist;
}
using namespace nanoflann;
void extract_patch_expand(DefaultMesh&i,arma::uvec&indices,DefaultMesh&o,int k)
{
    MeshKDTreeInterface<DefaultMesh> pts(i);
    KDTreeSingleIndexAdaptor<
            L2_Simple_Adaptor<float,MeshKDTreeInterface<DefaultMesh>>,
            MeshKDTreeInterface<DefaultMesh>,
            3,arma::uword>
            kdtree(3,pts,KDTreeSingleIndexAdaptorParams(5));
    kdtree.buildIndex();
    arma::uvec expand_indices(k);
    arma::fvec dists(k);
    std::vector<arma::uword> expanded_vec = arma::conv_to<std::vector<arma::uword>>::from(indices);
    float* points = (float*)i.points();
//    std::cerr<<"searching:"<<std::endl;
    expanded_vec.reserve(5000);
    for(int idx=0;idx<indices.n_rows;++idx)
    {
        kdtree.knnSearch(&points[3*indices(idx)],k,expand_indices.memptr(),dists.memptr());
        for(int n=0;n<k;++n)expanded_vec.push_back(expand_indices(n));
    }
    arma::uvec expanded_indices(expanded_vec.data(),expanded_vec.size(),false,true);
    arma::uvec  unique_indices = arma::unique(expanded_indices);
    extractMesh<DefaultMesh,DefaultMesh>(i,unique_indices,o);
}

void extract_patch_expand(DefaultMesh&i,arma::uvec&indices,DefaultMesh&o,float r)
{
    MeshKDTreeInterface<DefaultMesh> pts(i);
    KDTreeSingleIndexAdaptor<
            L2_Simple_Adaptor<float,MeshKDTreeInterface<DefaultMesh>>,
            MeshKDTreeInterface<DefaultMesh>,
            3,arma::uword>
            kdtree(3,pts,KDTreeSingleIndexAdaptorParams(5));
    kdtree.buildIndex();
    std::vector<std::pair<arma::uword,float>> expand_indices;
    std::vector<std::pair<arma::uword,float>>::iterator iiter;
    std::vector<arma::uword> expanded_vec = arma::conv_to<std::vector<arma::uword>>::from(indices);
    float* points = (float*)i.points();
    arma::fmat pmat(points,3,i.n_vertices(),false,true);
    arma::fmat ppmat = pmat.cols(indices);
    arma::fvec min = arma::min(ppmat,1);
//    std::cerr<<"searching:"<<std::endl;
    expanded_vec.reserve(5000);
    for(int idx=0;idx<indices.n_rows;++idx)
    {
        kdtree.radiusSearch(&points[3*indices(idx)],r,expand_indices,SearchParams());
        for(iiter=expand_indices.begin();iiter!=expand_indices.end();++iiter)
        {
            arma::fvec added = pmat.col(iiter->first);
            if(added(2)>min(2))
            {
                expanded_vec.push_back(iiter->first);
            }
        }
    }
    arma::uvec expanded_indices(expanded_vec.data(),expanded_vec.size(),false,true);
    arma::uvec  unique_indices = arma::unique(expanded_indices);
    extractMesh<DefaultMesh,DefaultMesh>(i,unique_indices,o);
}
