#include "extractpatchfeature.h"
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
