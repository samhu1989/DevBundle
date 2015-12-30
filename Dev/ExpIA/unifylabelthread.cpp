#include "unifylabelthread.h"
#include "extractpatchfeature.h"
#include "mbb.h"
#include <QTime>
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
    QTime timer;
    timer.restart();
    emit message(QString("Extracting Feature"),0);
    extract_patch_features();
    emit message(QString("Learning GMM model"),0);
    learn();
    emit message(QString("Assigning Patch Label"),0);
    assign();
    QString msg;
    int ms = timer.elapsed();
    int s = ms/1000;
    ms -= s*1000;
    int m = s/60;
    s -= m*60;
    int h = m/60;
    m -= h*60;
    msg = msg.sprintf("Time Used of Unify Label:%2u:%2u:%2u.%3u",h,m,s,ms);
    emit message(msg,0);
}

void UnifyLabelThread::extract_patch_features()
{
    input_patch_label_value_.clear();
    patch_features_.clear();
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
        ++lindex;
    }
    //reduce dimension
    int custom_dim = config_->getInt("Feature_dim");
    if(config_->has("Feature_dim")){
        if( feature_dim > custom_dim )
        {
            if(feature_base_.is_empty())
            {
                arma::mat feature_data;
                std::vector<arma::mat>::iterator fiter;
                for(fiter=patch_features_.begin();fiter!=patch_features_.end();++fiter)
                {
                    feature_data = arma::join_rows(feature_data,*fiter);
                }
                std::cerr<<"Feature dimension is "<<feature_data.n_rows<<std::endl;
                arma::vec mean = arma::mean(feature_data,1);
                //centralize data
                feature_data.each_col() -= mean;
                arma::mat U,V;
                arma::vec s;
                arma::svd(U,s,V,feature_data.t());
                emit message(tr("Init Feature Base with PCA"),0);
                std::cerr<<"Init Feature Base with PCA"<<std::endl;
                feature_base_ = V.cols(0,custom_dim-1);
            }
            std::vector<arma::mat>::iterator fiter;
            for(fiter=patch_features_.begin();fiter!=patch_features_.end();++fiter)
            {
//                std::cerr<<"n before reduce"<<(*fiter).n_cols<<std::endl;
                *fiter = ((*fiter).t()*feature_base_).t();
                if(custom_dim != (*fiter).n_rows)
                {
                    std::cerr<<"supposed to reduced to "<<custom_dim<<std::endl;
                    std::cerr<<"actually reduced to"<<(*fiter).n_rows<<std::endl;
                }
//                std::cerr<<"n after reduce"<<(*fiter).n_cols<<std::endl;
            }
        }
    }
    std::cerr<<"Feature dimension is reduced to"<<custom_dim<<std::endl;
}

void UnifyLabelThread::learn()
{
    size_t max_patch_num = 0;
    size_t max_patch_num_index = 0;
    size_t feature_dim = 2; //default to reduce dimension to 2
    std::vector<arma::mat>::iterator fiter;
    size_t index = 0;
    arma::mat data;
    for(fiter=patch_features_.begin();fiter!=patch_features_.end();++fiter)
    {
        if((*fiter).n_cols>max_patch_num)
        {
            max_patch_num=(*fiter).n_cols;
            max_patch_num_index = index;
        }
        data = arma::join_rows(data,*fiter);
        ++index;
    }
    feature_dim = data.n_rows;
    gmm_.reset(feature_dim,max_patch_num);
    std::cerr<<"choose frame "<<max_patch_num_index<<" with "<<max_patch_num<<"patches as clustering center"<<std::endl;
    gmm_.set_means(patch_features_[max_patch_num_index]);
    arma::mat dcovs(feature_dim,gmm_.n_gaus());
    dcovs.each_col() = arma::var(data,0,1);
    dcovs += std::numeric_limits<double>::epsilon();
    gmm_.set_dcovs(dcovs);
    gmm_.learn(data,max_patch_num,arma::maha_dist,arma::keep_existing,0,30,1e-10,true);
}

void UnifyLabelThread::assign()
{
    std::vector<arma::mat>::iterator fiter;
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
        const arma::mat& features,
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
    log_p.max(col_max_index,row_max_index);
    std::cerr<<"log_p:"<<std::endl;
    std::cerr<<log_p<<std::endl;
    while(!log_p.is_empty())
    {
        std::cerr<<gaus_indices(col_max_index)<<"->"<<feature_indices(row_max_index)<<":"<<log_p(col_max_index,row_max_index)<<std::endl;
        label_value( feature_indices(row_max_index) ) = 1 + arma::uword(gaus_indices(col_max_index));
        log_p.shed_row(col_max_index);
        gaus_indices.shed_row(col_max_index);
        log_p.shed_col(row_max_index);
        feature_indices.shed_row(row_max_index);
        if(log_p.is_empty())return;
        log_p.max(col_max_index,row_max_index);
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

