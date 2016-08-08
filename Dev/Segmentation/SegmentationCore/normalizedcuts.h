#ifndef NORMALIZEDCUTS_H
#define NORMALIZEDCUTS_H
#include "common.h"
#include <armadillo>
#include <QImage>
namespace Segmentation{
template<typename Mesh>
class NormalizedCuts
{
public:
    NormalizedCuts();
    bool configure(Config::Ptr);
    void cutW(std::shared_ptr<arma::sp_mat>w,arma::uvec&);
    void cutImage(const QImage&img, arma::uvec&);
    void cutMesh(typename MeshBundle<Mesh>::Ptr m,arma::uvec&);
    void cutGraph(typename MeshBundle<Mesh>::Ptr m,arma::uvec&);
    inline void getLabel(arma::uvec& label){label=label_;}
protected:
    void computeW_Image(const QImage& img);
    void computeW_Mesh(typename MeshBundle<Mesh>::Ptr m);
    void computeW_Graph(typename MeshBundle<Mesh>::Ptr m);
    void decompose();
    void clustering_GMM();
    void computeLabel_GMM();
    void bicut();
    void multicut();
    //for image
public:
    void createKernels(std::vector<arma::mat>&kernels);
protected:
    double scale0_;
    double scale1_;
    double scale2_;
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
private:
    std::shared_ptr<arma::sp_mat> W_;
    arma::sp_mat A_;
    arma::vec lambda_;
    arma::mat Y_;
    arma::vec th_;
    arma::uvec label_;
    double tol_;
    arma::uword k_;
    arma::uword max_N_;
    arma::gmm_diag gmm_;
};
}
#endif // NORMALIZEDCUTS_H
