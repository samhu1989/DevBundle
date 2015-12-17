#include "unifylabelthread.h"
#include "featurecore.h"
#include "mbb.h"
bool UnifyLabelThread::configure(Config::Ptr config)
{
    config_ = config;
    if(!config_->has("Color_Space"))return false;
    if("RGB"==config_->getString("Color_Space"))
    {
        if(!config_->has("RGB_r_bin_num"))return false;
        if(!config_->has("RGB_g_bin_num"))return false;
        if(!config_->has("RGB_b_bin_num"))return false;
    }else if("Lab"==config_->getString("Color_Space")){
        if(!config_->has("Lab_L_bin_num"))return false;
        if(!config_->has("Lab_a_bin_num"))return false;
        if(!config_->has("Lab_b_bin_num"))return false;
    }else return false;
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
    arma::mat dcovs(feature_dim,gmm_.n_gaus());
    dcovs.each_col() = arma::var(data,0,1);
    dcovs += std::numeric_limits<double>::epsilon();
    gmm_.set_dcovs(dcovs);
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
    arma::mat feature_mat = arma::conv_to<arma::mat>::from(features);
    uint64_t N_gaus = gmm_.n_gaus();
    uint64_t N_feature = feature_mat.n_cols;
    arma::mat log_p(N_gaus,N_feature);
    arma::uvec gaus_indices = arma::linspace<arma::uvec>(0,N_gaus-1,N_gaus);
    arma::uvec feature_indices = arma::linspace<arma::uvec>(0,N_feature-1,N_feature);
//    std::cerr<<"0"<<std::endl;
    label_value = arma::urowvec(N_feature,arma::fill::zeros);
//    std::cerr<<"1"<<std::endl;
    for( size_t gindex = 0 ; gindex < gmm_.n_gaus() ; ++gindex )
    {
        log_p.row(gindex) = gmm_.log_p(feature_mat,gindex);
    }
//    std::cerr<<"2"<<std::endl;
    uint64_t row_max_index;
    uint64_t col_max_index;
    uint64_t current_row = 0;
    log_p.max(col_max_index,row_max_index);
    current_row = col_max_index;
    while(!log_p.is_empty())
    {
        while(current_row >= log_p.n_rows ) -- current_row;
        log_p.row(current_row).max(row_max_index);
        log_p.col(row_max_index).max(col_max_index);
        if(current_row==col_max_index)
        {
            std::cerr<<gaus_indices(current_row)<<"->"<<feature_indices(row_max_index)<<std::endl;
            label_value( feature_indices(row_max_index) ) = 1 + arma::uword(gaus_indices(current_row));
//            std::cerr<<"a"<<std::endl;
            log_p.shed_row(current_row);
            gaus_indices.shed_row(current_row);
//            std::cerr<<"c"<<std::endl;
            log_p.shed_col(row_max_index);
//            std::cerr<<"d"<<std::endl;
            feature_indices.shed_row(row_max_index);
        }else{
            log_p.max(col_max_index,row_max_index);
            current_row = col_max_index;
        }
//        std::cerr<<current_row<<std::endl;
//        std::cerr<<"log_p.n_rows:"<<log_p.n_rows<<std::endl;
//        std::cerr<<"log_p.n_cols:"<<log_p.n_cols<<std::endl;
    }
//    std::cerr<<"3"<<std::endl;
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

