#include "featureview.h"
#include "ui_featureview.h"
#include "extractpatchfeature.h"
featureview::featureview(
        std::vector<MeshBundle<DefaultMesh>::Ptr>& inputs,
        std::vector<arma::uvec>& labels,
        arma::mat& base,
        QWidget *parent
        ):
    inputs_(inputs),
    labels_(labels),
    feature_base_(base),
    QWidget(parent),
    ui(new Ui::featureview)
{
    ui->setupUi(this);
    viewer_ = new FeatureViewerWidget();
    ui->gridLayout->addWidget(viewer_);
}

bool featureview::configure(Config::Ptr config)
{
    config_ = config;
}

void featureview::init(void)
{
    extract_patch_features();
    set_patch_features();
    viewer_->refresh();
}

void featureview::extract_patch_features()
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
        ++ lindex;
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
                emit message(tr("Reset Feature Base with PCA"),0);
                std::cerr<<"Reset Feature Base with PCA"<<std::endl;
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

void featureview::set_patch_features()
{
    arma::fmat features;
    std::vector<arma::mat>::iterator fiter;
    for(fiter=patch_features_.begin();fiter!=patch_features_.end();++fiter)
    {
        features = arma::join_rows(features,arma::conv_to<arma::fmat>::from(*fiter));
    }
    arma::Mat<uint8_t> colors(3,features.n_cols);
    QStringList strings;
    for(size_t fIdx=0;fIdx<input_patch_label_value_.size();++fIdx)
    {
        int pN = input_patch_label_value_[fIdx].size();
        QString str;
        for(size_t pIdx=0;pIdx<pN;++pIdx)
        {
            arma::uword label = input_patch_label_value_[fIdx](pIdx);

            str = str.sprintf("Frame %u Patch %u",fIdx,label);
            strings.push_back(str);
        }
    }
    viewer_->set_features(features);
    viewer_->set_feature_colors(colors);
    viewer_->set_feature_string(strings);
}

featureview::~featureview()
{
    viewer_->deleteLater();
    delete ui;
}
