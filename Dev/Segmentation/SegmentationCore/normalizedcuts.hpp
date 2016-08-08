#include "normalizedcuts.h"
namespace Segmentation{
template<typename Mesh>
NormalizedCuts<Mesh>::NormalizedCuts()
{
    k_ = 30;
    tol_ = 0.01;
    scale0_ = 1.0;
    scale1_ = 2.0;
    scale2_ = 2.0;
    kernel_size_ = 7;
    max_N_ = 10;
}
template<typename Mesh>
bool NormalizedCuts<Mesh>::configure(Config::Ptr config)
{
    ;
}
template<typename Mesh>
void NormalizedCuts<Mesh>::cutW(std::shared_ptr<arma::sp_mat>w,arma::uvec&label)
{
    W_ = w;
    decompose();
    clustering_GMM();
    computeLabel_GMM();
    getLabel(label);
}
template<typename Mesh>
void NormalizedCuts<Mesh>::cutImage(const QImage&img,arma::uvec&label)
{
    computeW_Image(img);
    decompose();
    clustering_GMM();
//    computeLabel_GMM();
    bicut();
    getLabel(label);
}
template<typename Mesh>
void NormalizedCuts<Mesh>::cutMesh(typename MeshBundle<Mesh>::Ptr m,arma::uvec&label)
{
    computeW_Mesh(m);
    decompose();
    clustering_GMM();
    computeLabel_GMM();
    getLabel(label);
}
template<typename Mesh>
void NormalizedCuts<Mesh>::cutGraph(typename MeshBundle<Mesh>::Ptr m,arma::uvec&label)
{
    computeW_Graph(m);
    decompose();
    clustering_GMM();
    computeLabel_GMM();
    getLabel(label);
}

double cut(arma::sp_mat& A, arma::vec& v, double threshold)
{
    double cut=0.0;
    for(size_t i=0;i<v.size();i++)
    {
        bool lower=( v(i) < threshold );
        for(size_t j=0;j<v.size();j++)
        {
            if(lower && v(j)>=threshold )
                cut+=A(i,j);
        }
    }
    return cut;
}

double assoc(arma::sp_mat& A, arma::vec& v, double threshold,bool a)
{
    double cut=0.0;
    for(size_t i=0;i<v.size();i++)
    {
        bool lower;
        if(a)
            lower=(v(i)<threshold);
        else
            lower=(v(i)>=threshold);
        if(lower)
        {
            for(size_t j=0;j<v.size();j++)
            {
                cut+=A(i,j);
            }
        }
    }
    return cut;
}

template<typename Mesh>
void NormalizedCuts<Mesh>::multicut()
{
//    arma::vec v = Y_.col(0);
//    label_ = arma::uvec(v.size(),arma::fill::zeros);
//    #pragma omp parallel for
//    for(size_t i=0;i<v.size();i++)
//    {
//        if(v(i)<best_threshold)
//        {
//            label_(i)=1;
//        }
//        else
//        {
//            label_(i)=2;
//        }
//    }
}

template<typename Mesh>
void NormalizedCuts<Mesh>::bicut()
{
    std::cerr<<"find best threshold"<<std::endl;
    double best_threshold=0.0;
//    double lowest_cost=std::numeric_limits<double>::max();
    arma::vec v = Y_.col(0);
//    float min = v.min();
//    float max = v.max();
//    float step = ( max - min ) / 10.0;
//    for(float i = min ;i < max ;i+=step)
//    {
//        std::cerr<<"i:"<<i<<std::endl;
//        double threshold=i;
//        double cutAB = cut(*W_,v,threshold);
//        std::cerr<<"cutAB:"<<cutAB<<std::endl;
//        double assocAV = assoc(*W_,v,threshold,true);
//        std::cerr<<"assocAV:"<<assocAV<<std::endl;
//        double assocBV = assoc(*W_,v,threshold,false);
//        std::cerr<<"assocBV:"<<assocBV<<std::endl;
//        double cost = cutAB/(1.0+assocAV) + cutAB/(1.0+assocBV);
//        std::cerr<<"cost:"<<cost<<std::endl;
//        if(cost<lowest_cost)
//        {
//            lowest_cost=cost;
//            best_threshold=threshold;
//        }
//    }
    std::cerr<<"best_threshold:"<<std::endl;
    std::cerr<<best_threshold<<std::endl;
    //create Mask
    label_ = arma::uvec(v.size(),arma::fill::zeros);
    #pragma omp parallel for
    for(size_t i=0;i<v.size();i++)
    {
        if(v(i)<best_threshold)
        {
            label_(i)=1;
        }
        else
        {
            label_(i)=2;
        }
    }
}

template<typename Mesh>
void NormalizedCuts<Mesh>::getGaussianKernel2D(double sigma1,double sigma2,double angle, int size, arma::mat& kernel)
{
    assert(size%2==1);
    arma::vec center(2);
    center.fill(size/2+1.0);
    arma::mat rot = getRotationMatrix2D(center,angle,1.0);
    arma::mat sigma(3,2);
    sigma(0,0)=1.0/sigma1;
    sigma(0,1)=0;
    sigma(1,0)=0;
    sigma(1,1)=1.0/sigma2;
    sigma(2,0)=0.0;
    sigma(2,1)=0.0;

    arma::mat invsigma=rot*sigma;
    double factor = 1.0/(2.0*M_PI*sqrt(sigma1*sigma2));
    arma::vec mu(2);
    mu(0)=size/2+1.0;
    mu(1)=size/2+1.0;

    kernel = arma::mat(size,size,arma::fill::zeros);
    arma::vec x(2);
    #pragma omp parallel for
    for(arma::uword c=0;c<size;c++)
    {
        double * col =(double*)kernel.colptr(c);
        for(arma::uword r=0;r<size;r++)
        {
            x(0)=r;
            x(1)=c;
            arma::vec tmp = invsigma*(x-mu);
            double m=arma::dot((x-mu),(tmp));
            col[r]=factor*exp(-0.5*m);
        }
    }
}

template<typename Mesh>
void NormalizedCuts<Mesh>::createKernels(std::vector<arma::mat>&kernels)
{
    kernels.resize(2);
    #pragma omp for
    for(uint32_t idx=0;idx<2;++idx)
    {
        getGaussianKernel2D(1.0+0.5*idx,0.5+0.5*idx,0.0,kernel_size_,kernels[idx]);
    }
//    #pragma omp for
//    for(uint32_t idx=4;idx<kernels.size();++idx)
//    {
//        getGaussianKernel2D(2.0,0.3,40.0*(idx-4),kernel_size_,kernels[idx]);
//    }
}

void blur(const arma::mat& imgmat,size_t width,size_t height,const arma::mat& kernel,arma::mat& blured)
{
    for(uint32_t idx=0;idx<imgmat.n_rows;++idx)
    {
        arma::mat reshaped = arma::reshape(imgmat.row(idx),width,height);
        arma::mat B = arma::conv2(reshaped,kernel,"same");
        blured.col(idx) = arma::vectorise(B);
    }
}

template<typename Mesh>
void NormalizedCuts<Mesh>::blurImage(const QImage& img,arma::mat& blured)
{
    arma::Mat<uint8_t> imgMat((uint8_t*)img.bits(),4,img.byteCount()/4,false,true);
    arma::Mat<uint8_t> tmp  = imgMat.rows(0,2);
    arma::fmat ftmp;
    ColorArray::RGB2Lab(tmp,ftmp);
    arma::mat imgDMat = arma::conv_to<arma::mat>::from(ftmp);
    std::vector<arma::mat> kernels;
    std::cerr<<"createKernels"<<std::endl;
    createKernels(kernels);
    blured = arma::mat(imgDMat.n_cols,imgDMat.n_rows*kernels.size(),arma::fill::zeros);
    double* m_ptr = (double*)blured.memptr();
    uint64_t offset = imgDMat.n_cols*imgDMat.n_rows;
    std::cerr<<"bluring"<<std::endl;
    #pragma omp parallel for
    for(uint32_t idx=0;idx<kernels.size();++idx)
    {
        arma::mat bluredbykernel((double*)(m_ptr+idx*offset),imgDMat.n_cols,imgDMat.n_rows,false,true);
        blur(imgDMat,img.width(),img.height(),kernels[idx],bluredbykernel);
    }
    arma::inplace_trans(blured);
}

template<typename Mesh>
void NormalizedCuts<Mesh>::computeW_Image(const QImage& img)
{
    QImage img8888 = img.convertToFormat(QImage::Format_RGBA8888);
    arma::Mat<uint8_t> imgMat((uint8_t*)img8888.bits(),4,img8888.byteCount()/4,false,true);
    arma::Mat<uint8_t> tmp  = imgMat.rows(0,2);
    arma::fmat ftmp;
    ColorArray::RGB2Lab(tmp,ftmp);
    arma::mat imgDMat = arma::conv_to<arma::mat>::from(ftmp);
    arma::mat bluredMat;
    std::cerr<<"blurImage"<<std::endl;
    blurImage(img8888,bluredMat);
    assert(bluredMat.is_finite());
//    bluredMat.save("./debug/label/bluredMat.arma",arma::raw_ascii);
    std::cerr<<"End blurImage"<<std::endl;
    size_t N = img.width()*img.height();
    W_.reset(new arma::sp_mat(N,N));
    *W_ = arma::speye(N,N);
//    imgDMat.save("./debug/label/imgDMat.arma",arma::raw_ascii);
    for(int r=0;r<img8888.height();r++)
    {
        std::cerr<<"r:"<<r<<std::endl;
        for(int c=0;c<img8888.width();c++)
        {
            uint32_t wi = r*img8888.width()+c;
            for(int i=std::max(0,r-2);i<std::min(img8888.height(),r+2);i++)
            {
                for(int j=std::max(0,c-2);j<std::min(img8888.width(),c+2);j++)
                {
                    uint32_t wj = i*img8888.width()+j;
                    if( (*W_)(wi,wj) != 0.0 )continue;
                    double affinity = distanceAffinity(r,c,i,j,scale0_);
                    affinity += vecAffinity<arma::vec>(imgDMat.col(wi),imgDMat.col(wj),scale1_);
//                    affinity += vecAffinity<arma::vec>(bluredMat.col(wi),bluredMat.col(wj),scale2_);
                    (*W_)(wi,wj) = std::exp(affinity);
                    (*W_)(wj,wi) = (*W_)(wi,wj);
                }
            }
        }
    }
    std::cerr<<"End W"<<std::endl;
//    arma::mat(*W_).save("./debug/label/W_.arma",arma::raw_ascii);
}
template<typename Mesh>
void NormalizedCuts<Mesh>::computeW_Mesh(typename MeshBundle<Mesh>::Ptr m)
{
    ;
}
template<typename Mesh>
void NormalizedCuts<Mesh>::computeW_Graph(typename MeshBundle<Mesh>::Ptr m)
{
    ;
}
template<typename Mesh>
void NormalizedCuts<Mesh>::decompose()
{
    std::cerr<<"decompose:"<<std::endl;
    arma::vec D = arma::vectorise(arma::mat(arma::sum(*W_)));
    arma::vec sqrt_D = arma::sqrt(D);
    arma::vec inv_sqrt_D = 1.0 / sqrt_D;
    arma::sp_mat Dmat = arma::speye<arma::sp_mat>(W_->n_rows,W_->n_cols);
    arma::sp_mat inv_sqrt_Dmat = arma::speye<arma::sp_mat>(W_->n_rows,W_->n_cols);
    Dmat.diag() = D;
    inv_sqrt_Dmat.diag() = inv_sqrt_D;
    A_ = inv_sqrt_Dmat*( Dmat - (*W_) )*inv_sqrt_Dmat;
    A_.diag() += 1000.0*std::numeric_limits<float>::epsilon();
//    A_.save("./debug/label/A.arma");
    bool success = false;
    double tol = 50.0;
    for(uint32_t try_n=0;try_n<7;++try_n)
    {
        std::cerr<<"try_n:"<<try_n<<std::endl;
        std::cerr<<"try tol:"<<tol<<std::endl;
        if(arma::eigs_sym(lambda_,Y_,A_,(k_+1),"sm",tol)){
            success = true;
            tol /= 10.0;
        }else break;
    }
    if(!success)std::cerr<<"Failed on decomposition, Please relax the tol"<<std::endl;//failed
//    lambda_.save("./debug/label/lambda_.arma",arma::raw_ascii);
//    Y_.save("./debug/label/Y_.arma",arma::raw_ascii);
    Y_.shed_col(0);
    Y_.each_col() %= inv_sqrt_D;
}
template<typename Mesh>
void NormalizedCuts<Mesh>::clustering_GMM()
{
    arma::mat Yt = Y_.t();
    arma::mat last_means;
    arma::mat last_covs;
    arma::mat last_hefts;
    arma::mat means = arma::mean(Yt,1);
    arma::mat covs(k_,1);
    std::cerr<<Yt.n_rows<<","<<Yt.n_cols<<std::endl;
    covs.col(0) = arma::var(Yt,0,1);
    arma::rowvec hefts(1,1,arma::fill::ones);
    double last_log_p = std::numeric_limits<double>::max();
    double log_p = 0.5*tol_*last_log_p;
    uint32_t iter = 0;
    while( ( ( log_p - last_log_p ) / last_log_p ) > tol_ || iter < 2 )
    {
        last_means = means;
        last_covs = covs;
        last_hefts = hefts;
        last_log_p = log_p;
        //increase k and update means covs hefts
        arma::rowvec maxcov = arma::mean(covs);
        arma::uword insert;
        maxcov.max(insert);
        arma::vec insert_mean(last_means.n_rows,arma::fill::randn);
        insert_mean *= 10.0*std::numeric_limits<float>::epsilon();
        insert_mean +=  last_means.col(insert);
        means = arma::join_rows(last_means,insert_mean);
        covs = arma::join_rows(last_covs,last_covs.col(insert));
        hefts = arma::rowvec(means.n_cols);
        hefts.fill(1.0/double(means.n_cols));
        //
//        gmm_.reset(means.n_rows,means.n_cols);
        gmm_.set_params(means,covs,hefts);
        gmm_.learn(Yt,means.n_cols,arma::eucl_dist,arma::keep_existing,0,10,1e-10,true);
        if(gmm_.n_gaus()>max_N_){
            std::cerr<<"NormalizedCuts<Mesh>::clustering():reach max_N:"<<max_N_<<std::endl;
            break;
        }
        means = gmm_.means;
        covs = gmm_.dcovs;
        hefts = gmm_.hefts;
        log_p = gmm_.avg_log_p(Yt);
        ++iter ;
        std::cerr<<"N:"<<means.n_cols<<std::endl;
        std::cerr<<"max_cov:"<<log_p<<std::endl;
        std::cerr<<"last_max_cov:"<<last_log_p<<std::endl;
        std::cerr<<"criteria:"<<( log_p - last_log_p ) / last_log_p  <<std::endl;
    }
}
template<typename Mesh>
void NormalizedCuts<Mesh>::computeLabel_GMM()
{
    label_ = (gmm_.assign(Y_.t(),arma::eucl_dist)).t();
    label_ += 1;
}
}
