#include "normalizedcuts.h"
namespace Segmentation{
template<typename Mesh>
NormalizedCuts<Mesh>::NormalizedCuts()
{
    k_ = 30;
    tol_ = 0.001;
    scale0_ = 1.0;
    scale1_ = 1.0;
    scale2_ = 1.0;
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
    clustering();
    computeLabel();
    getLabel(label);
}
template<typename Mesh>
void NormalizedCuts<Mesh>::cutImage(const QImage&img,arma::uvec&label)
{
    std::cerr<<"computeW_Image(img)"<<std::endl;
    computeW_Image(img);
    std::cerr<<"decompose"<<std::endl;
    decompose();
    std::cerr<<"clustering"<<std::endl;
    clustering();
    std::cerr<<"computeLabel"<<std::endl;
    computeLabel();
    std::cerr<<"getLabel(label)"<<std::endl;
    getLabel(label);
}
template<typename Mesh>
void NormalizedCuts<Mesh>::cutMesh(typename MeshBundle<Mesh>::Ptr m,arma::uvec&label)
{
    computeW_Mesh(m);
    decompose();
    clustering();
    computeLabel();
    getLabel(label);
}
template<typename Mesh>
void NormalizedCuts<Mesh>::cutGraph(typename MeshBundle<Mesh>::Ptr m,arma::uvec&label)
{
    computeW_Graph(m);
    decompose();
    clustering();
    computeLabel();
    getLabel(label);
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
    kernels.resize(8);
    #pragma omp for
    for(uint32_t idx=0;idx<4;++idx)
    {
        getGaussianKernel2D(0.5+0.5*idx,0.5+0.5*idx,0.0,kernel_size_,kernels[idx]);
    }
    #pragma omp for
    for(uint32_t idx=4;idx<kernels.size();++idx)
    {
        getGaussianKernel2D(2.0,0.3,40.0*(idx-4),kernel_size_,kernels[idx]);
    }
}

void blur(const arma::mat& imgmat,size_t width,size_t height,const arma::mat& kernel,arma::mat& blured)
{
    for(uint32_t idx=0;idx<imgmat.n_rows;++idx)
    {
        arma::mat reshaped = arma::reshape(imgmat.row(idx),height,width);
        arma::mat B = arma::conv2(reshaped,kernel,"same");
        blured.col(idx) = arma::vectorise(B);
    }
}

template<typename Mesh>
void NormalizedCuts<Mesh>::blurImage(const QImage& img,arma::mat& blured)
{
    arma::Mat<uint8_t> imgMat((uint8_t*)img.bits(),3,img.width()*img.height(),false,true);
    arma::mat imgDMat = arma::conv_to<arma::mat>::from(imgMat);
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
    QImage img888 = img.convertToFormat(QImage::Format_RGB888);
    arma::Mat<uint8_t> imgMat((uint8_t*)img888.bits(),3,img888.width()*img888.height(),false,true);
    arma::mat bluredMat;
    std::cerr<<"blurImage"<<std::endl;
    blurImage(img888,bluredMat);
    std::cerr<<"End blurImage"<<std::endl;
    size_t N = img.width()*img.height();
    W_.reset(new arma::sp_mat(N,N));
    for(int r=0;r<img888.height();r++)
    {
        for(int c=0;c<img888.width();c++)
        {
            uint32_t wi = r*img.width()+c;
            for(int i=std::max(0,r-5);i<std::min(img888.height(),r+5);i++)
            {
                for(int j=std::max(0,c-5);j<std::min(img888.width(),c+5);j++)
                {
                    uint32_t wj = i*img888.width()+j;
                    double affinity = distanceAffinity(r,c,i,j,scale0_);
                    affinity += vecAffinity<arma::Col<uint8_t>>(imgMat.col(wi),imgMat.col(wj),scale1_);
                    affinity += vecAffinity<arma::vec>(bluredMat.col(wi),bluredMat.col(wj),scale2_);
                    (*W_)(wi,wj) = std::exp(affinity);
                    (*W_)(wj,wi) = (*W_)(wi,wj);
                }
            }
        }
    }
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
    arma::vec D = arma::vectorise(arma::mat(arma::sum(*W_)));
    arma::vec sqrt_D = arma::sqrt(D);
    arma::vec inv_sqrt_D = 1.0 / sqrt_D;
    arma::sp_mat Dmat = arma::speye<arma::sp_mat>(W_->n_rows,W_->n_cols);
    arma::sp_mat inv_sqrt_Dmat = arma::speye<arma::sp_mat>(W_->n_rows,W_->n_cols);
    Dmat.diag() = D;
    inv_sqrt_Dmat.diag() = inv_sqrt_D;
    A_ = inv_sqrt_Dmat*( Dmat - (*W_) )*inv_sqrt_Dmat;
    assert(arma::eigs_sym(lambda_,Y_,A_,(k_+1),"sm"));
    Y_.shed_col(0);
    Y_.each_col() %= sqrt_D;
}
template<typename Mesh>
void NormalizedCuts<Mesh>::clustering()
{
    double last_max_cov = std::numeric_limits<double>::max();
    double max_cov = 0.5*tol_*last_max_cov;
    arma::mat Yt = Y_.t();
    arma::mat last_means;
    arma::mat last_covs;
    arma::mat last_hefts;
    arma::mat means = arma::mean(Yt,1);
    arma::mat covs(k_,1);
    std::cerr<<Yt.n_rows<<","<<Yt.n_cols<<std::endl;
    covs.col(0) = arma::var(Yt,0,1);
    arma::rowvec hefts(1,1,arma::fill::ones);
    while( ( max_cov / last_max_cov ) < tol_ )
    {
        std::cerr<<"max_cov:"<<max_cov<<std::endl;
        std::cerr<<"last_max_cov:"<<last_max_cov<<std::endl;
        last_means = means;
        last_covs = covs;
        last_hefts = hefts;
        last_max_cov = max_cov;
        //increase k and update means covs hefts
        arma::rowvec maxcov = arma::max(covs);
        arma::uword insert;
        maxcov.max(insert);
        arma::vec insert_mean(last_means.n_rows,arma::fill::randn);
        insert_mean *= 2.0*std::numeric_limits<float>::epsilon();
        insert_mean +=  last_means.col(insert);
        means = arma::join_rows(last_means,insert_mean);
        covs = arma::join_rows(last_covs,last_covs.col(insert));
        hefts = arma::rowvec(means.n_cols);
        hefts.fill(1.0/double(means.n_cols));
        //
//        gmm_.reset(means.n_rows,means.n_cols);
        gmm_.set_params(means,covs,hefts);
        gmm_.learn(Yt,means.n_cols,arma::eucl_dist,arma::keep_existing,30,0,1e-10,true);
        if(gmm_.n_gaus()>max_N_){
            std::cerr<<"NormalizedCuts<Mesh>::clustering():reach max_N:"<<max_N_<<std::endl;
            break;
        }
        means = gmm_.means;
        covs = gmm_.dcovs;
        hefts = gmm_.hefts;
        max_cov = arma::max(arma::max(covs));
    }
}
template<typename Mesh>
void NormalizedCuts<Mesh>::computeLabel()
{
    label_ = (gmm_.assign(Y_.t(),arma::eucl_dist)).t();
    label_ += 1;
}
}
