#include "unifylabelthread.h"
#include "featurecore.h"
#include "mbb.h"
bool UnifyLabelThread::configure(Config::Ptr config)
{
    config_ = config;
    if(!config_->has("RGB_r_bin_num"))return false;
    if(!config_->has("RGB_g_bin_num"))return false;
    if(!config_->has("RGB_b_bin_num"))return false;
    return true;
}

void UnifyLabelThread::run()
{
    emit message(QString("Extracting Feature"),0);
    extract_patch_features();
    emit message(QString("Learning GMM model"),0);
    learn();
    emit message(QString("Assigning Patch Label"),0);
    assign();
    emit message(QString("Done Unify Label"),0);
}

void UnifyLabelThread::extract_patch_features()
{
    input_patch_label_value_.clear();
    patch_features_.clear();
    MeshBundle<DefaultMesh>::PtrList::iterator iter;
    size_t lindex = 0;
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
                arma::fvec feature;
                extract_patch_feature(extracted_mesh,feature);
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
        ++lindex;
    }
}

void UnifyLabelThread::extract_patch_feature(DefaultMesh&mesh,arma::fvec&feature)
{
    arma::fmat points((float*)mesh.points(),3,mesh.n_vertices(),false,true);
    arma::fmat box;
    get3DMBB(points,2,box);
    arma::fvec lwh(3);//length width height
    lwh(0) = arma::norm( box.col(0) - box.col(1) );
    lwh(1) = arma::norm( box.col(0) - box.col(3) );
    lwh(2) = arma::norm( box.col(0) - box.col(4) );
    if( lwh(0) < lwh(1) ) lwh.swap_rows(0,1);
    arma::fvec color_hist;
    Feature::ColorHistogramRGB<DefaultMesh>
            color_hist_extractor(
                config_->getInt("RGB_r_bin_num"),
                config_->getInt("RGB_g_bin_num"),
                config_->getInt("RGB_b_bin_num")
                );
    color_hist_extractor.extract(mesh,color_hist);
    if(config_->has("Feature_height_w"))
    {
        lwh.tail(1) *= config_->getFloat("Feature_height_w");
    }
    if(config_->has("Feature_length_width_w"))
    {
        lwh.head(2) *= config_->getFloat("Feature_length_width_w");
    }
    if(config_->has("Feature_color_hist_w"))
    {
        color_hist *= config_->getFloat("Feature_color_hist_w");
    }
    feature = arma::join_cols(lwh,color_hist);
//    feature = color_hist;
}

void UnifyLabelThread::learn()
{
    size_t max_patch_num = 0;
    size_t max_patch_num_index = 0;
    size_t feature_dim = 0;
    std::vector<arma::fmat>::iterator fiter;
    size_t index = 0;
    arma::mat data;
    for(fiter=patch_features_.begin();fiter!=patch_features_.end();++fiter)
    {
        if(0==feature_dim)
        {
            feature_dim = (*fiter).n_rows;
        }
        if((*fiter).n_cols>max_patch_num)
        {
            max_patch_num=(*fiter).n_cols;
            max_patch_num_index = index;
        }
        data = arma::join_rows(data,arma::conv_to<arma::mat>::from(*fiter));
        ++index;
    }
    gmm_.reset(feature_dim,max_patch_num);
    gmm_.set_means(arma::conv_to<arma::mat>::from(patch_features_[max_patch_num_index]));
    gmm_.learn(data,max_patch_num,arma::eucl_dist,arma::keep_existing,10,0,1e-10,true);
}

void UnifyLabelThread::assign()
{
    std::vector<arma::fmat>::iterator fiter;
    size_t index = 0;
    for(fiter=patch_features_.begin();fiter!=patch_features_.end();++fiter)
    {
        std::cerr<<"assigning F"<<index<<std::endl;
        arma::urowvec output_patch_label_value;
        assign(*fiter,output_patch_label_value);
        alter_label(input_patch_label_value_[index],output_patch_label_value,labels_[index]);
        MeshBundle<DefaultMesh>::Ptr iptr = inputs_[index];
        iptr->custom_color_.fromlabel(labels_[index]);
        ++index;
    }
}

void UnifyLabelThread::assign(
        const arma::fmat& features,
        arma::urowvec& label_value
        )
{
    label_value = arma::urowvec(features.n_cols,arma::fill::zeros);
    arma::rowvec label_value_p(features.n_cols);
    label_value_p.fill( std::numeric_limits<double>::lowest() );
    arma::mat feature_mat = arma::conv_to<arma::mat>::from(features);
    for( size_t gindex = 0 ; gindex < gmm_.n_gaus() ; ++gindex )
    {
        arma::rowvec log_p = gmm_.log_p(feature_mat,gindex);
        arma::uvec log_p_indices = arma::linspace<arma::uvec>(0,log_p.size()-1,log_p.size());
        arma::uword max_p_index;
        bool assigned = false;
        while(!assigned)
        {
            log_p.max(max_p_index);
            if( log_p(max_p_index) > label_value_p( log_p_indices(max_p_index) ) )
            {
                label_value( log_p_indices(max_p_index) ) = gindex + 1;
                label_value_p( log_p_indices(max_p_index) ) = log_p(max_p_index);
                assigned = true;
            }else{
                log_p.shed_col(max_p_index);
                log_p_indices.shed_row(max_p_index);
            }
            if(log_p.is_empty())break;
        }
//        std::cerr<<"gindex:"<<std::endl;
//        std::cerr<<gindex<<std::endl;
//        if(!assigned)std::cerr<<"not assinged"<<std::endl;
//        else std::cerr<<"assigned to "<<log_p_indices(max_p_index)<<std::endl;
//        std::cerr<<"log_p:"<<std::endl;
//        std::cerr<<log_p<<std::endl;
    }
}

void UnifyLabelThread::alter_label(
        const arma::urowvec& in_label_value,
        const arma::urowvec& out_label_value,
        arma::uvec &label
        )
{
    assert(in_label_value.size()==out_label_value.size());
//    std::cerr<<"in_label_value:"<<std::endl;
//    std::cerr<<in_label_value<<std::endl;
//    std::cerr<<"in_label_value.size():"<<std::endl;
//    std::cerr<<in_label_value.size()<<std::endl;
//    std::cerr<<"out_label_value:"<<std::endl;
//    std::cerr<<out_label_value<<std::endl;
//    std::cerr<<"out_label_value.size():"<<std::endl;
//    std::cerr<<out_label_value.size()<<std::endl;
    arma::uvec altered_label = label;
    for(size_t index=0 ; index < in_label_value.size() ; ++index )
    {
        arma::uvec indices = arma::find( label == in_label_value[index] );
        if(!indices.is_empty()){
            altered_label(indices).fill( out_label_value[index] );
        }else{
            std::cerr<<"indices empty according to in_label_value"<<std::endl;
        }
    }
    label = altered_label;
}

