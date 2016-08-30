#include "normalizedcuts.h"
#include <QTime>
#include <QColor>
namespace Segmentation{
template<typename Mesh>
NormalizedCuts<Mesh>::NormalizedCuts()
{
    type_ = N;
    k_ = 40;
    clustering_k_ = 5;
    tol_ = 0.01;
    d_scale_ = 1.0;
    c_scale_ = 1.0;
    f_scale_ = 1.0;
    convex_scale_ = 10.0*std::numeric_limits<float>::epsilon();
    kernel_size_ = 7;
    max_N_ = 10;
    eps_ = 0.0;
}
template<typename Mesh>
bool NormalizedCuts<Mesh>::configure(Config::Ptr config)
{
    if(config->has("NCut_Type"))
    {
        if(config->getString("NCut_Type")=="Min"||config->getString("NCut_Type")=="min")type_=M;
        if(config->getString("NCut_Type")=="GSP"||config->getString("NCut_Type")=="gsp")type_=G;
    }
    if(config->has("NCut_Clustering"))
    {
        if(config->getString("NCut_Clustering")=="GMM")clustering_type_=GMM;
        if(config->getString("NCut_Clustering")=="Bisection")clustering_type_=Bisection;
        if(config->getString("NCut_Clustering")=="Kmean")clustering_type_=Kmean;
    }
    if(config->has("NCut_Clustering_k"))
    {
        clustering_k_ = config->getInt("NCut_Clustering_k");
    }
    if(config->has("NCut_eps")){
        eps_ = config->getDouble("NCut_eps");
    }
    if(config->has("NCut_k")){
        k_ = config->getInt("NCut_k");
    }
    if(config->has("NCut_GMM_tol"))
    {
        tol_ = config->getDouble("NCut_GMM_tol");
    }
    if(config->has("NCut_GMM_Max_n"))
    {
        max_N_ = config->getInt("NCut_GMM_Max_n");
    }
    if(config->has("NCut_Dist_Scale"))
    {
        d_scale_ = config->getDouble("NCut_Dist_Scale");
    }
    if(config->has("NCut_Color_Scale"))
    {
        c_scale_ = config->getDouble("NCut_Color_Scale");
    }
    if(config->has("NCut_Feature_Scale"))
    {
        f_scale_ = config->getDouble("NCut_Feature_Scale");
    }
    if(config->has("NCut_Convexity_Scale"))
    {
        convex_scale_ = config->getDouble("NCut_Convexity_Scale");
    }
    return true;
}
template<typename Mesh>
void NormalizedCuts<Mesh>::cutW(std::shared_ptr<arma::sp_mat>w,arma::uvec&label)
{
    W_ = w;
    decompose();
    clustering();
    getLabel(label);
}
template<typename Mesh>
void NormalizedCuts<Mesh>::cutImage(const QImage&img,arma::uvec&label)
{
    computeW_Image(img);
    decompose();
    clustering();
    getLabel(label);
}
template<typename Mesh>
void NormalizedCuts<Mesh>::cutMesh(typename MeshBundle<Mesh>::Ptr m,arma::uvec&label)
{
    computeW_Mesh(m);
    decompose();
    clustering();
    getLabel(label);
}
template<typename Mesh>
void NormalizedCuts<Mesh>::cutGraph(typename MeshBundle<Mesh>::Ptr m,arma::uvec&label)
{
    computeW_Graph(m);
    decompose();
    clustering();
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
    throw std::logic_error("void NormalizedCuts<Mesh>::multicut(): not implemented");
}

template<typename Mesh>
void NormalizedCuts<Mesh>::bicut()
{
    std::cerr<<"Bisection"<<std::endl;
    double best_threshold=0.0;
    arma::vec v = Y_.col(0);
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
    kernels.resize(6);
    #pragma omp parallel for
    for(uint32_t idx=0;idx<4;++idx)
    {
        getGaussianKernel2D(1.0+0.5*idx,1.0+0.5*idx,0.0,kernel_size_,kernels[idx]);
    }
    #pragma omp parallel for
    for(uint32_t idx=4;idx<6;++idx)
    {
        getGaussianKernel2D(2.0,1.5,idx*90.0,kernel_size_,kernels[idx]);
    }
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
    arma::mat imgDMat = arma::conv_to<arma::mat>::from(tmp);
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
    arma::mat imgDMat = arma::conv_to<arma::mat>::from(imgMat.rows(0,2));
    arma::mat bluredMat;
    std::cerr<<"blurImage"<<std::endl;
    blurImage(img8888,bluredMat);
    assert(bluredMat.is_finite());
    std::cerr<<"End blurImage"<<std::endl;
    size_t N = img.width()*img.height();
    W_.reset(new arma::sp_mat(N,N));
    *W_ = arma::speye(N,N);
    for(int r=0;r<img8888.height();r++)
    {
//        std::cerr<<"r:"<<r<<std::endl;
        for(int c=0;c<img8888.width();c++)
        {
            uint32_t wi = r*img8888.width()+c;
            for(int i=std::max(0,r-3);i<std::min(img8888.height(),r+3);i++)
            {
                for(int j=std::max(0,c-3);j<std::min(img8888.width(),c+3);j++)
                {
                    uint32_t wj = i*img8888.width()+j;
                    if( (*W_)(wi,wj) != 0.0 )continue;
                    double affinity = distanceAffinity(r,c,i,j,d_scale_);
                    affinity += vecAffinity<arma::vec>(imgDMat.col(wi),imgDMat.col(wj),c_scale_);
//                    affinity += vecAffinity<arma::vec>(bluredMat.col(wi),bluredMat.col(wj),f_scale_);
                    (*W_)(wi,wj) = std::exp(affinity);
                    (*W_)(wj,wi) = (*W_)(wi,wj);
                }
            }
        }
    }
    std::cerr<<"End W"<<std::endl;
}

template<typename Mesh>
void NormalizedCuts<Mesh>::computeW_Mesh(typename MeshBundle<Mesh>::Ptr m)
{
    ;
}

template<typename Mesh>
void NormalizedCuts<Mesh>::computeW_Graph(typename MeshBundle<Mesh>::Ptr m)
{
    MeshBundle<Mesh>& mesh = *m;
    VoxelGraph<Mesh>& graph = mesh.graph_;
    size_t N = graph.size();
    W_.reset(new arma::sp_mat(N,N));
    *W_ = arma::speye(N,N);
    for(arma::Mat<uint16_t>::iterator niter=graph.voxel_neighbors.begin();niter!=graph.voxel_neighbors.end();   )
    {
        uint16_t wi = *niter;
        ++niter;
        uint16_t wj = *niter;
        ++niter;
        double affinity = vecAffinity<arma::fvec>(
                    graph.voxel_centers.col(wi),
                    graph.voxel_centers.col(wj),
                    d_scale_
                    );
        affinity = std::exp(affinity);
        affinity += convexity<arma::fvec>(
                    graph.voxel_centers.col(wi),
                    graph.voxel_normals.col(wi),
                    graph.voxel_centers.col(wj),
                    graph.voxel_normals.col(wj),
                    convex_scale_
                    );
        (*W_)(wi,wj) = 0.5*affinity;
        (*W_)(wj,wi) = (*W_)(wi,wj);
    }
}

template<typename Mesh>
void NormalizedCuts<Mesh>::debug_convexity(typename MeshBundle<Mesh>::Ptr m)
{
    MeshBundle<Mesh>& mesh = *m;
    VoxelGraph<Mesh>& graph = mesh.graph_;
    size_t N = graph.size();
    graph.voxel_edge_colors = arma::Mat<uint8_t>(4,graph.voxel_neighbors.n_cols,arma::fill::zeros);
    arma::vec conv(graph.voxel_edge_colors.n_cols,arma::fill::zeros);
    arma::vec::iterator conv_iter = conv.begin();
    for(arma::Mat<uint16_t>::iterator niter=graph.voxel_neighbors.begin();niter!=graph.voxel_neighbors.end();   )
    {
        uint16_t wi = *niter;
        ++niter;
        uint16_t wj = *niter;
        ++niter;
        double affinity = convexity<arma::fvec>(
                    graph.voxel_centers.col(wi),
                    graph.voxel_normals.col(wi),
                    graph.voxel_centers.col(wj),
                    graph.voxel_normals.col(wj),
                    convex_scale_
                    );
        *conv_iter = affinity;
        ++conv_iter;
    };
    //to do convert convexity to color
    double max = arma::max(conv);
    double min = arma::min(conv);
    std::cerr<<"Convexity["<<min<<","<<max<<"]"<<std::endl;
    conv -= min;
    conv /= ( max - min );
    assert(conv.is_finite());
    uint32_t* ptr = (uint32_t*)graph.voxel_edge_colors.memptr();
    ColorArray::colorfromValue(ptr,graph.voxel_edge_colors.n_cols,conv);
}

template<typename Mesh>
void NormalizedCuts<Mesh>::debug_color(typename MeshBundle<Mesh>::Ptr m)
{
    MeshBundle<Mesh>& mesh = *m;
    VoxelGraph<Mesh>& graph = mesh.graph_;
    size_t N = graph.size();
    graph.voxel_edge_colors = arma::Mat<uint8_t>(4,graph.voxel_neighbors.n_cols,arma::fill::zeros);
    arma::fmat Lab;
    ColorArray::RGB2Lab(graph.voxel_colors,Lab);
    arma::fmat ab = Lab.rows(0,1);
    arma::vec wv(graph.voxel_edge_colors.n_cols,arma::fill::zeros);
    arma::vec::iterator wv_iter = wv.begin();
    for(arma::Mat<uint16_t>::iterator niter=graph.voxel_neighbors.begin();niter!=graph.voxel_neighbors.end();   )
    {
        uint16_t wi = *niter;
        ++niter;
        uint16_t wj = *niter;
        ++niter;
        double affinity = vecAffinity<arma::fvec>(
                    ab.col(wi),
                    ab.col(wj),
                    c_scale_
                    );
        *wv_iter = std::exp(affinity);
        ++wv_iter;
    };
    //to do convert convexity to color
    double max = arma::max(wv);
    double min = arma::min(wv);
    std::cerr<<"Color["<<min<<","<<max<<"]"<<std::endl;
    wv -= min;
    wv /= ( max - min );
    assert(wv.is_finite());
    uint32_t* ptr = (uint32_t*)graph.voxel_edge_colors.memptr();
    ColorArray::colorfromValue(ptr,graph.voxel_edge_colors.n_cols,wv);
}
template<typename Mesh>
void NormalizedCuts<Mesh>::debug_dist(typename MeshBundle<Mesh>::Ptr m)
{
    MeshBundle<Mesh>& mesh = *m;
    VoxelGraph<Mesh>& graph = mesh.graph_;
    size_t N = graph.size();
    graph.voxel_edge_colors = arma::Mat<uint8_t>(4,graph.voxel_neighbors.n_cols,arma::fill::zeros);
    arma::vec wv(graph.voxel_edge_colors.n_cols,arma::fill::zeros);
    arma::vec::iterator wv_iter = wv.begin();
    for(arma::Mat<uint16_t>::iterator niter=graph.voxel_neighbors.begin();niter!=graph.voxel_neighbors.end();   )
    {
        uint16_t wi = *niter;
        ++niter;
        uint16_t wj = *niter;
        ++niter;
        double affinity = vecAffinity<arma::fvec>(
                    graph.voxel_centers.col(wi),
                    graph.voxel_centers.col(wj),
                    d_scale_
                    );
        *wv_iter = std::exp(affinity);
        ++wv_iter;
    };
    //to do convert convexity to color
    double max = arma::max(wv);
    double min = arma::min(wv);
    std::cerr<<"Dist["<<min<<","<<max<<"]"<<std::endl;
    wv -= min;
    wv /= ( max - min );
    assert(wv.is_finite());
    uint32_t* ptr = (uint32_t*)graph.voxel_edge_colors.memptr();
    ColorArray::colorfromValue(ptr,graph.voxel_edge_colors.n_cols,wv);
}
template<typename Mesh>
void NormalizedCuts<Mesh>::debug_W(typename MeshBundle<Mesh>::Ptr m)
{
    MeshBundle<Mesh>& mesh = *m;
    VoxelGraph<Mesh>& graph = mesh.graph_;
    size_t N = graph.size();
    graph.voxel_edge_colors = arma::Mat<uint8_t>(4,graph.voxel_neighbors.n_cols,arma::fill::zeros);
    arma::vec wv(graph.voxel_edge_colors.n_cols,arma::fill::zeros);
    arma::vec::iterator wv_iter = wv.begin();
    for(arma::Mat<uint16_t>::iterator niter=graph.voxel_neighbors.begin();niter!=graph.voxel_neighbors.end();   )
    {
        uint16_t wi = *niter;
        ++niter;
        uint16_t wj = *niter;
        ++niter;
        double affinity = vecAffinity<arma::fvec>(
                    graph.voxel_centers.col(wi),
                    graph.voxel_centers.col(wj),
                    d_scale_
                    );
        affinity = std::exp(affinity);
        affinity += convexity<arma::fvec>(
                    graph.voxel_centers.col(wi),
                    graph.voxel_normals.col(wi),
                    graph.voxel_centers.col(wj),
                    graph.voxel_normals.col(wj),
                    convex_scale_
                    );
        *wv_iter = affinity;
        ++wv_iter;
    };
    //to do convert convexity to color
    double max = arma::max(wv);
    double min = arma::min(wv);
    std::cerr<<"W=["<<min<<","<<max<<"]"<<std::endl;
    wv -= min;
    wv /= ( max - min );
    assert(wv.is_finite());
    uint32_t* ptr = (uint32_t*)graph.voxel_edge_colors.memptr();
    ColorArray::colorfromValue(ptr,graph.voxel_edge_colors.n_cols,wv);
}

template<typename Mesh>
void NormalizedCuts<Mesh>::decomposeNormarlized()
{
    std::cerr<<"NormalizedCut:"<<std::endl;
    arma::vec D = arma::vectorise(arma::mat(arma::sum(*W_)));
    arma::uvec zeroIndex = arma::find(D <= 0);
    if(!zeroIndex.empty())
    {
        zeroIndex.print("zeroIndex:");
    }
    arma::vec sqrt_D = arma::sqrt(D);
    arma::vec inv_sqrt_D = 1.0 / sqrt_D;
    arma::sp_mat Dmat = arma::speye<arma::sp_mat>(W_->n_rows,W_->n_cols);
    arma::sp_mat inv_sqrt_Dmat = arma::speye<arma::sp_mat>(W_->n_rows,W_->n_cols);
    Dmat.diag() = D;
    inv_sqrt_Dmat.diag() = inv_sqrt_D;
    A_ = inv_sqrt_Dmat*( Dmat - (*W_) )*inv_sqrt_Dmat;
//    QTime time;
//    arma::mat(A_).save("./debug/label/A_"+time.currentTime().toString("hh_mm_ss").toStdString()+".arma",arma::raw_ascii);
    bool success = false;
    double stol = 50.0;
    double etol = std::numeric_limits<double>::epsilon();
    success = arma_custom::eigs_sym(lambda_,Y_,A_,(k_+1),"sm",stol,etol);
    if(!success)std::cerr<<"Failed on decomposition, Please relax the tol"<<std::endl;//failed
//    lambda_.print("lambda_:");
    std::cerr<<"eps:"<<eps_<<std::endl;
    arma::uvec index = arma::find( lambda_ <= eps_ );
//    index.print("index:");
    if(!index.empty())Y_.shed_cols(index(0),index(index.size()-1));
    Y_.each_col() %= inv_sqrt_D;
}

template<typename Mesh>
void NormalizedCuts<Mesh>::decomposeMin()
{
    std::cerr<<"MinCut:"<<std::endl;
    arma::vec D = arma::vectorise(arma::mat(arma::sum(*W_)));
    arma::sp_mat Dmat = arma::speye<arma::sp_mat>(W_->n_rows,W_->n_cols);
    Dmat.diag() = D;
    A_ = Dmat - (*W_);
//    double eps = std::min(std::abs(arma::min(arma::min(A_))),std::abs(arma::max(arma::max(A_))));
//    double eps = std::numeric_limits<float>::epsilon();
//    QTime time;
//    arma::mat(A_).save("./debug/label/A_"+time.currentTime().toString("hh_mm_ss").toStdString()+".arma",arma::raw_ascii);
    bool success = false;
    double stol = 50.0;
    double etol = std::numeric_limits<double>::epsilon();
    success = arma_custom::eigs_sym(lambda_,Y_,A_,(k_+1),"sm",stol,etol);
    if(!success)std::cerr<<"Failed on decomposition, Please relax the tol"<<std::endl;//failed
//    lambda_.print("lambda_:");
    arma::uvec index = arma::find( lambda_ <= eps_ );
    std::cerr<<"eps:"<<eps_<<std::endl;
//    index.print("index:");
    if(!index.empty())Y_.shed_cols(index(0),index(index.size()-1));
}

template<typename Mesh>
void NormalizedCuts<Mesh>::decomposeGSP()
{
    std::cerr<<"GSP:"<<std::endl;
    arma::vec D = arma::vectorise(arma::mat(arma::sum(*W_)));
    arma::sp_mat Dmat = arma::speye<arma::sp_mat>(W_->n_rows,W_->n_cols);
    Dmat.diag() = D;
    A_ = Dmat - (*W_);
    bool success = false;
    double stol = 50.0;
    double etol = std::numeric_limits<double>::epsilon();
    success = arma_custom::eigs_sym(lambda_,Y_,A_,(k_+1),"sm",stol,etol);
    if(!success)std::cerr<<"Failed on decomposition, Please relax the tol"<<std::endl;//failed
//    lambda_.print("lambda_:");
    arma::uvec index = arma::find( lambda_ <= 0.0 );
//    index.print("index:");
    if(!index.empty()){
        lambda_.shed_rows(index(0),index(index.size()-1));
        Y_.shed_cols(index(0),index(index.size()-1));
    }
    arma::vec sqrt_lambda_ = arma::sqrt(lambda_);
    Y_.each_row() /= sqrt_lambda_.t();
}
template<typename Mesh>
void NormalizedCuts<Mesh>::clustering_GMM()
{
    std::cerr<<"GMM"<<std::endl;
    arma::mat Yt = Y_.t();
    arma::uword k = Yt.n_rows;
    arma::mat last_means;
    arma::mat last_covs;
    arma::mat last_hefts;
    arma::mat means = arma::mean(Yt,1);
    arma::mat covs(k,1);
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
        gmm_.set_params(means,covs,hefts);
        gmm_.learn(Yt,means.n_cols,arma::eucl_dist,arma::keep_existing,0,50,1e-10,false);
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
template<typename Mesh>
void NormalizedCuts<Mesh>::clustering_Kmean()
{
    std::cerr<<"Kmean"<<std::endl;
    std::cerr<<"Clustering_k_:"<<clustering_k_<<std::endl;
    arma::vec kratio(1);
    std::mt19937 engine;  // Mersenne twister random number engine
    std::uniform_real_distribution<double> distr(0.0, 1.0);
    kratio.imbue( [&]() { return distr(engine); } );
    arma::uword k = kratio(0) * double(clustering_k_) + clustering_k_;
    arma::mat Yt = Y_.t();
    gmm_.learn(Yt,k,arma::eucl_dist,arma::random_spread,50,0,1e-10,true);
}
template<typename Mesh>
void NormalizedCuts<Mesh>::computeLabel_Kmean()
{
    label_ = (gmm_.assign(Y_.t(),arma::eucl_dist)).t();
    label_ += 1;
}
}
