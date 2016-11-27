#ifndef NORMALIZEDCUTS_H
#define NORMALIZEDCUTS_H
#include "common.h"
#include <armadillo>
#include <QImage>
#include <random>
namespace Segmentation{
template<typename Mesh>
class NormalizedCuts
{
public:
    typedef enum{
        N,//normalized cut
        M, //min cut
        G  //Global Point Signature
    }TYPE;
    typedef enum{
        GMM,
        Kmean,
        Bisection
    }CLUSTERING;
public:
    NormalizedCuts();
    bool configure(Config::Ptr);
    void cutW(std::shared_ptr<arma::sp_mat>w,arma::uvec&);
    void cutImage(const QImage&img, arma::uvec&);
    void cutMesh(typename MeshBundle<Mesh>::Ptr m,arma::uvec&);
    void cutGraph(typename MeshBundle<Mesh>::Ptr m,arma::uvec&);

    inline void getLabel(arma::uvec& label){label=label_;}
    inline void setK(const arma::uword& k){k_=k;}
    inline void setType(const typename NormalizedCuts<Mesh>::TYPE t){type_=t;}
    inline const arma::mat& getY()const{return Y_;}

    void computeW_Image(const QImage& img);
    void computeW_Mesh(typename MeshBundle<Mesh>::Ptr m);
    void computeW_Graph(typename MeshBundle<Mesh>::Ptr m);

    void decompose()
    {
        switch(type_)
        {
        case N:decomposeNormarlized();break;
        case M:decomposeMin();break;
        case G:decomposeGPS();break;
        default:decomposeNormarlized();
        }
    }
protected:
    void decomposeNormarlized();
    void decomposeMin();
    void decomposeGPS();
    void clustering()
    {
        switch(clustering_type_)
        {
        case GMM:
            clustering_GMM();
            computeLabel_GMM();
            break;
        case Bisection:bicut();break;
        case Kmean:
            clustering_Kmean();
            computeLabel_Kmean();
            break;
        default:bicut();
        }
    }
    void clustering_GMM();
    void computeLabel_GMM();
    void clustering_Kmean();
    void computeLabel_Kmean();
    void bicut();
    void multicut();
    //for image
public:
    void createKernels(std::vector<arma::mat>&kernels);
protected:
    double d_scale_;
    double c_scale_;
    double f_scale_;
    double convex_scale_;
    uint32_t kernel_size_;
    inline double distanceAffinity(double x1, double y1, double x2, double y2,double scale)
    {
        return -((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2))/scale;
    }
    template<typename vec>
    inline double vecAffinity(const vec& x,const vec& y,double scale)
    {
        vec dif = x-y;
        return - arma::dot(dif,dif)/scale;
    }
    void getGaussianKernel2D(double sigma1, double sigma2, double angle, int size, arma::mat &kernel);
    void blurImage(const QImage&,arma::mat&);
    //for graph and mesh
public:
    void debug_convexity(typename MeshBundle<Mesh>::Ptr m);
    void debug_color(typename MeshBundle<Mesh>::Ptr m);
    void debug_dist(typename MeshBundle<Mesh>::Ptr m);
    void debug_W(typename MeshBundle<Mesh>::Ptr m);
protected:
    template<typename vec>
    inline double convexity(const vec& p0,const vec&n0,const vec&p1,const vec& n1,double eps)
    {
        eps = eps > 0?100.0*eps:std::numeric_limits<double>::epsilon();
        vec dir = arma::normalise(p1-p0);
        double dif = std::acos(arma::dot(n0,dir)) - std::acos(arma::dot(n1,dir));
        if(dif<-eps)return 0.0;
        else if(dif>eps)return 1.0;
        else return 0.5*dif/eps+0.5;
    }
protected:
    TYPE type_;
    CLUSTERING  clustering_type_;
    std::shared_ptr<arma::sp_mat> W_;
    arma::sp_mat A_;
    arma::vec lambda_;
    arma::mat Y_;
    arma::vec th_;
    arma::uvec label_;
    double tol_;
    double eps_;
    arma::uword k_;
    arma::uword clustering_min_k_;
    arma::uword clustering_max_k_;
    arma::uword max_N_;
    arma::gmm_diag gmm_;
    std::mt19937 rand_engine_;  // Mersenne twister random number engine
    std::uniform_int_distribution<arma::uword> distr_;
};
}
#endif // NORMALIZEDCUTS_H
