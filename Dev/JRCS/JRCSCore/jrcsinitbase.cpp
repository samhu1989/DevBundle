#include "jrcsinitbase.h"
#include "featurecore.h"
JRCSInitBase::JRCSInitBase()
{
    ;
}


void JRCSInitBase::init_alpha_with_label(
        MatPtrLst& alpha,
        const MatPtrLst& vv,
        const MatPtrLst& vn,
        const CMatPtrLst& vc,
        const LCMatPtrLst& vlc,
        const LMatPtrLst& vl,
        bool verbose
        )
{
    ;
}

bool JRCSInitBase::configure(Config::Ptr config)
{
    config_ = config;
    return true;
}

void JRCSInitBase::extract_patch_features()
{
    input_patch_label_value_.clear();
    patch_features_.clear();
    size_t lindex = 0;
    int feature_dim;
    for( lindex = 0 ; lindex < vv_.size() ; ++lindex )
    {
        arma::uword label_value = 1;
        arma::uword label_max = arma::max((*vl_[lindex]));
        while( label_value <= label_max )
        {
            arma::uvec extracted_label = arma::find( label_value == (*vl_[lindex]) );
            DefaultMesh extracted_mesh;
            if(!extracted_label.is_empty())
            {
                extractMesh<DefaultMesh>(*vv_[lindex],*vn_[lindex],*vc_[lindex],extracted_label,extracted_mesh);
                arma::vec feature;
                extract_patch_feature(extracted_mesh,feature,config_);
                feature_dim = feature.size();
                if( input_patch_label_value_.size() <= lindex )
                {
                    input_patch_label_value_.emplace_back(1);
                    input_patch_label_value_.back()[0] = label_value;
                    patch_features_.emplace_back(feature.size(),1);
                    patch_features_.back().col(0) = feature;
                }else{
                    input_patch_label_value_.back().insert_cols(0,1);
                    input_patch_label_value_.back()[0] = label_value;
                    patch_features_.back().insert_cols(0,feature);
                }
            }
            ++ label_value;
        }
    }
}

void JRCSInitBase::learn()
{
    ;
}

void JRCSInitBase::assign()
{
    ;
}
