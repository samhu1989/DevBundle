#include "jrcsinitbase.h"
#include "featurecore.h"
#include <armadillo>
namespace JRCS
{
JRCSInitBase::JRCSInitBase()
{
    ;
}

bool JRCSInitBase::init_with_label(
        const int k,
        const MatPtrLst& vv,
        const MatPtrLst& vn,
        const CMatPtrLst& vc,
        const LCMatPtrLst& vlc,
        const LMatPtrLst& vl,
        int verbose
        )
{
    k_ = k;
    vv_ = vv;
    vn_ = vn;
    vc_ = vc;
    vl_ = vl;
    vlc_ = vlc;
    verbose_ = verbose;

    if(verbose_)std::cerr<<"Init by "<<name()<<std::endl;

    if(verbose_>0)std::cerr<<"extract_patch_features()"<<std::endl;
    extract_patch_features();
    if(verbose_>0)std::cerr<<"pca()"<<std::endl;
    pca();
    if(verbose_>0)std::cerr<<"learn()"<<std::endl;
    learn();
    if(verbose_>0)std::cerr<<"assign()"<<std::endl;
    assign();
    if(verbose_>0)std::cerr<<"generate_prob()"<<std::endl;
    generate_prob();
    if(verbose_>0)std::cerr<<"generate_alpha()"<<std::endl;
    generate_alpha();
    return true;
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
                if(verbose_>0)
                {
                    std::cerr<<"start extract_patch_feature"<<std::endl;
                    std::cerr<<"label_value:"<<label_value<<std::endl;
                    std::cerr<<"extracted_mesh.n_vertices():"<<extracted_mesh.n_vertices()<<std::endl;
                }
                extract_patch_feature(extracted_mesh,feature,config_);
                if(verbose_>0)
                {
                    std::cerr<<"done extract_patch_feature"<<std::endl;
                }
                assert(feature.is_finite());
//                feature_dim = feature.size();
                if( input_patch_label_value_.size() <= lindex )
                {
                    input_patch_label_value_.emplace_back(1);
                    input_patch_label_value_.back()[0] = label_value;
                    patch_features_.emplace_back(feature.size(),1);
                    patch_features_.back().col(0) = feature;
                    patch_sizes_.emplace_back(1);
                    patch_sizes_.back().row(0) = extracted_mesh.n_vertices();
                }else{
                    input_patch_label_value_.back().insert_cols(0,1);
                    input_patch_label_value_.back()[0] = label_value;
                    patch_features_.back().insert_cols(0,feature);
                    patch_sizes_.back().insert_rows(0,1);
                    patch_sizes_.back().row(0) = extracted_mesh.n_vertices();
                }
            }
            ++ label_value;
        }
    }
}

void JRCSInitBase::pca()
{
    int custom_dim = config_->getInt("Feature_dim");
    int feature_dim = patch_features_[0].n_rows;
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
//                std::cerr<<"Feature dimension is "<<feature_data.n_rows<<std::endl;
                arma::vec mean = arma::mean(feature_data,1);
                //centralize data
                arma::mat centered_feature = feature_data.each_col() - mean;
                arma::mat V;
                arma::vec s;
                arma::mat A = (centered_feature*centered_feature.t());
                arma::eig_sym(s,V,A);
                if(verbose_>0){
                    std::cerr<<"using PCA"<<std::endl;
                    std::cerr<<A;
                    std::cerr<<std::endl;
                    std::cerr<<"s:"<<std::endl;
                    std::cerr<<s<<std::endl;
                }
//                emit message(tr("Init Feature Base with PCA"),0);
//                std::cerr<<"Init Feature Base with PCA"<<std::endl;
                arma::uvec sort_index = arma::sort_index(s,"descend");
                feature_base_ = arma::join_rows(mean,V.cols(sort_index.head(custom_dim)));
            }
            std::vector<arma::mat>::iterator fiter;
            arma::vec mean = feature_base_.col(0);
            arma::mat proj = feature_base_.cols(1,feature_base_.n_cols-1);
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
//    std::cerr<<"Feature dimension is reduced to"<<custom_dim<<std::endl;
}

bool JRCSInitBase::check_centers()
{
    int custom_dim = config_->getInt("Feature_dim");
    if(feature_centers_.n_rows!=custom_dim)return false;
    if(feature_centers_.n_cols<2)return false;
    return true;
}

void JRCSInitBase::learn()
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
    if(!check_centers()){
        if(verbose_>0)std::cerr<<"choose frame "<<max_patch_num_index<<" with "<<max_patch_num<<"patches as clustering center"<<std::endl;
        feature_centers_ = patch_features_[max_patch_num_index];
    }
    gmm_.reset(feature_centers_.n_rows,feature_centers_.n_cols);
    gmm_.set_means(feature_centers_);
    arma::mat dcovs(feature_dim,gmm_.n_gaus());
    dcovs.each_col() = arma::var(data,0,1);
    dcovs += std::numeric_limits<double>::epsilon();
    gmm_.set_dcovs(dcovs);
    gmm_.learn(data,max_patch_num,arma::eucl_dist,arma::keep_existing,20,0,1e-10,true);
    feature_centers_ = gmm_.means;
}

void JRCSInitBase::assign()
{
    std::vector<arma::mat>::iterator fiter;
    size_t index = 0;
    patch_prob_.resize(patch_features_.size());
    for(fiter=patch_features_.begin();fiter!=patch_features_.end();++fiter)
    {
        assign(*fiter,patch_prob_[index]);
        ++index;
    }
}

void JRCSInitBase::assign(const arma::mat& f,arma::fmat& p)
{
    uint32_t N_feature = f.n_cols;
    p = arma::fmat(gmm_.n_gaus(),f.n_cols,arma::fill::zeros);
    if(verbose_>0)std::cerr<<"p("<<p.n_rows<<","<<p.n_cols<<")"<<std::endl;
    for(size_t i = 0 ; i < gmm_.n_gaus() ; ++i )
    {
        arma::rowvec pi(N_feature);
        #pragma omp parallel for
        for(size_t findex = 0; findex < N_feature ; ++findex)
        {
            pi(findex) = - arma::norm( gmm_.means.col(i) - f.col(findex) );
        }
        p.row(i) = arma::conv_to<arma::frowvec>::from( arma::exp( pi ) );
    }
    if(verbose_>0)std::cerr<<"done assigning"<<std::endl;
}

void JRCSInitBase::generate_alpha()
{
    alpha_.resize(vv_.size());
    size_t index;
    DMatPtrLst::iterator iter;
    std::vector<arma::urowvec>::iterator liter = input_patch_label_value_.begin();
    arma::uword r_k = k_;
    arma::uvec obj_size(prob_.size());
    arma::fvec prob = prob_;
    //determine the gaussian numbers for each object
    for(int oi = 0 ; oi < prob_.size() ; ++ oi )
    {
        obj_size(oi) = arma::uword(float(k_)*float(prob(oi)));
        obj_size(oi) = std::max(arma::uword(9),obj_size(oi));
        obj_size(oi) = std::min(r_k,obj_size(oi));
        r_k -= obj_size(oi);
    }
    //calculate the column index by gaussian numbers
    std::vector<arma::uvec> cols;
    arma::uword sum = 0;
    for(int oi = 0 ; oi < obj_size.size() ; ++oi )
    {
        cols.emplace_back(obj_size[oi],arma::fill::zeros);
        cols.back() = arma::linspace<arma::uvec>(sum,sum+obj_size[oi]-1,obj_size[oi]);
        sum += obj_size[oi];
    }
    //fill alpha
    index = 0;
    for( iter = alpha_.begin() ; iter != alpha_.end() ; ++iter )
    {
        iter->reset(new arma::mat(vv_[index]->n_cols,k_,arma::fill::zeros));
        arma::mat& alpha = **iter;
        arma::fmat& prob = patch_prob_[index];
        arma::uvec& label = *vl_[index];
        alpha.fill(1.0/float(k_));
        #pragma omp parallel for
        for(int c=0;c<prob.n_cols;++c)
        {
            arma::uvec rows = arma::find(label==input_patch_label_value_[index](c));
            for(int r=0;r<prob.n_rows;++r)
            {
                alpha(rows,cols[r]).fill(prob(r,c));
            }
        }
        ++index;
    }
}

void JRCSInitBase::getAlpha(DMatPtrLst& alpha)
{
    alpha = alpha_;
}

void JRCSInitBase::generate_prob()
{
    std::vector<arma::uvec>::iterator iter;
    size_t index = 0;
    for(iter=patch_sizes_.begin();iter!=patch_sizes_.end();++iter)
    {
        arma::fmat& p = patch_prob_[index];
        arma::fvec size = arma::conv_to<arma::fvec>::from(*iter);
        if(prob_.n_rows!=p.n_rows)
        {
            prob_ = arma::fvec(p.n_rows,arma::fill::zeros);
        }
        prob_ += (p*size);
        ++index;
    }
    float sum = arma::accu(prob_);
    prob_ /= sum;
}

void JRCSInitBase::getObjProb(arma::fvec& obj_prob)
{
    obj_prob = prob_;
}

}
