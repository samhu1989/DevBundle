#include "featureview.h"
#include "ui_featureview.h"
#include "extractpatchfeature.h"
featureview::featureview(
        std::vector<MeshBundle<DefaultMesh>::Ptr>& inputs,
        std::vector<arma::uvec>& labels,
        arma::mat& base,
        arma::mat& center,
        QWidget *parent
        ):
    inputs_(inputs),
    labels_(labels),
    feature_base_(base),
    feature_center_(center),
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
    return true;
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
                arma::mat centered_feature = feature_data.each_col() - mean;
                arma::mat U,V;
                arma::vec s;
                arma::svd(U,s,V,centered_feature.t());
                emit message(tr("Init Feature Base with PCA"),0);
                std::cerr<<"Init Feature Base with PCA"<<std::endl;
                feature_base_ = V.cols(0,custom_dim-1);
                feature_base_ = arma::join_rows(mean,feature_base_);
                emit message(tr("Reset Feature Base with PCA"),0);
                std::cerr<<"Reset Feature Base with PCA"<<std::endl;
            }
            std::vector<arma::mat>::iterator fiter;
            arma::vec mean = feature_base_.col(0);
            arma::mat proj = feature_base_.cols(1,feature_base_.n_cols-1);
//            std::cerr<<mean<<std::endl;
            for(fiter=patch_features_.begin();fiter!=patch_features_.end();++fiter)
            {
//                std::cerr<<"n before reduce"<<(*fiter).n_cols<<std::endl;
                *fiter = (( (*fiter).each_col() - mean ).t()*proj).t();
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
    if(!feature_center_.is_empty())
    {
        features = arma::join_rows(arma::conv_to<arma::fmat>::from(feature_center_),features);
    }
    QStringList strings;
    arma::Mat<uint8_t> colors;
    colors = arma::Mat<uint8_t>( 3 , features.n_cols );
    std::cerr<<feature_center_.n_cols<<feature_center_.n_cols<<std::endl;
    for(size_t cIdx=0 ; cIdx < feature_center_.n_cols ; ++ cIdx )
    {
        QString str;
        colors(0,cIdx) = 0;
        colors(1,cIdx) = 0;
        colors(2,cIdx) = 0;
        str = str.sprintf("Object %u Centroid",cIdx+1);
        strings.push_back(str);
    }

    int index;
    size_t idx = feature_center_.n_cols;
    for(size_t fIdx=0;fIdx<input_patch_label_value_.size();++fIdx)
    {
        int pN = input_patch_label_value_[fIdx].size();
        QString str;
        for(size_t pIdx=0 ; pIdx < pN; ++pIdx)
        {
            arma::uword label = input_patch_label_value_[fIdx](pIdx);
            if(label==0)index=0;
            else{
                std::srand(label);
                index = std::rand()%( ColorArray::DefaultColorNum_ - 1 );
                index += 1;
            }
            colors(0,idx) = ColorArray::DefaultColor[index].rgba.r;
            colors(1,idx) = ColorArray::DefaultColor[index].rgba.g;
            colors(2,idx) = ColorArray::DefaultColor[index].rgba.b;
            str = str.sprintf("Frame %u Patch %u",fIdx,label);
            strings.push_back(str);
            ++idx;
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
