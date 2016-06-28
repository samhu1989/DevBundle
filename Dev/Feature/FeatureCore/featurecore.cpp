#include "colorhistogram.hpp"
#include "blockbasedfeature.hpp"
#include "featurecore.h"
#include "pointnormal.h"
void extract_patch_feature(DefaultMesh&mesh,arma::vec&feature,Config::Ptr config_)
{
//    std::cerr<<"0"<<std::endl;
    assert(mesh.n_vertices()>0);
    arma::vec block_feature;//length width height
//    std::vector<double> color_feature(125);
    {
        Feature::BlockBasedFeature<DefaultMesh> block_feature_extractor;
        block_feature_extractor.extract(mesh,block_feature);
    }
    assert(block_feature.is_finite());
    //
//    std::cerr<<"1"<<std::endl;
    if("RGB"==config_->getString("Color_Space"))
    {
        Feature::ColorHistogramRGB<DefaultMesh>
                color_hist_extractor(
                    config_->getInt("RGB_r_bin_num"),
                    config_->getInt("RGB_g_bin_num"),
                    config_->getInt("RGB_b_bin_num")
                    );
        color_hist_extractor.extract(mesh,feature);
    }else if("Lab"==config_->getString("Color_Space")){
        Feature::ColorHistogramLab<DefaultMesh>
                color_hist_extractor(
                    config_->getInt("Lab_L_bin_num"),
                    config_->getInt("Lab_a_bin_num"),
                    config_->getInt("Lab_b_bin_num")
                    );
        color_hist_extractor.extract(mesh,feature);
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
//    std::cerr<<"3"<<std::endl;
//    feature = arma::vec(color_feature);
    assert(feature.is_finite());
    if(config_->has("Feature_color_hist_w"))
    {
        feature *= config_->getFloat("Feature_color_hist_w");
    }
    feature = arma::join_cols(feature,block_feature);
//    feature = color_hist;
}
