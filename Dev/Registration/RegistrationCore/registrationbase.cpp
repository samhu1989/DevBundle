#include "RegistrationBase.h"
#include "common.h"
#include "nanoflann.hpp"
namespace Registration{
void RegistrationBase::alignAroundZ(const arma::fmat &x0,const arma::fmat &x1, arma::fmat &bestR)
{
    arma::fmat::fixed<3,3> R;
    arma::fmat y;
    float N = 100.0;
    double minError = std::numeric_limits<double>::max();
    for(int k=0;k<N;++k)
    {
        getRotation(0,0,0,0,0,k*2.0*M_PI/N,R);
        y = R*x1;
        double e = closestError(x0,y);
        if( minError >  e )
        {
            minError = e;
            bestR = R;
        }
    }
//    std::cerr<<"minError:"<<minError<<std::endl;
}
using namespace nanoflann;
double RegistrationBase::closestError(const arma::fmat&x0,const arma::fmat&x1)
{
    ArmaKDTreeInterface<arma::fmat> points(x0);
    KDTreeSingleIndexAdaptor<
            L2_Simple_Adaptor<float,ArmaKDTreeInterface<arma::fmat>>,
            ArmaKDTreeInterface<arma::fmat>,
            3,arma::uword>
            kdtree(3,points,KDTreeSingleIndexAdaptorParams(9));
    kdtree.buildIndex();
    double error = 0;
    float* pts = (float*)x1.memptr();
    arma::uvec indices(1);
    arma::fvec dists(1);
    for(size_t index=0;index<x1.n_cols;++index)
    {
        kdtree.knnSearch(&pts[3*index],1,indices.memptr(),dists.memptr());
        error += std::sqrt(dists(0));
    }
    error /= x1.n_cols;
    return error;
}
void RegistrationBase::alignAroundZ(
        const arma::fmat& x0,
        const arma::fmat& n0,
        const arma::fmat& x1,
        const arma::fmat& n1,
        arma::fmat& bestR
        )
{
    arma::fmat::fixed<3,3> R;
    arma::fmat y;
    arma::fmat yn;
    float N = 100.0;
    double minError = std::numeric_limits<double>::max();
    for(int k=0;k<N;++k)
    {
        getRotation(0,0,0,0,0,k*2.0*M_PI/N,R);
        y = R*x1;
        yn = R*n1;
        double e = closestError(x0,n0,y,yn);
        if( minError >  e )
        {
            minError = e;
            bestR = R;
        }
    }
//    std::cerr<<"minError:"<<minError<<std::endl;
}
double RegistrationBase::closestError(
        const arma::fmat& x,
        const arma::fmat& xn,
        const arma::fmat& y,
        const arma::fmat& yn
        )
{
    if(x.n_cols!=xn.n_cols||x.n_rows!=xn.n_rows||y.n_cols!=yn.n_cols||y.n_rows!=yn.n_rows)
    {
        throw std::logic_error("x and xn or y and yn are not the same size");
    }

    ArmaKDTreeInterface<arma::fmat> points(x);
    KDTreeSingleIndexAdaptor<
            L2_Simple_Adaptor<float,ArmaKDTreeInterface<arma::fmat>>,
            ArmaKDTreeInterface<arma::fmat>,
            3,arma::uword>
            kdtree(3,points,KDTreeSingleIndexAdaptorParams(9));
    kdtree.buildIndex();
    double error = 0;
    float* pts = (float*)y.memptr();
    arma::uvec indices(1);
    arma::fvec dists(1);
    for(size_t index=0;index<y.n_cols;++index)
    {
        kdtree.knnSearch(&pts[3*index],1,indices.memptr(),dists.memptr());
        float dcos = 0.0;
        if( yn.col(index).is_finite() && xn.col(indices(0)).is_finite() )dcos = arma::dot(yn.col(index),xn.col(indices(0)));
        error += ( 1.0 - dcos )*std::sqrt(dists(0));
    }
    error /= y.n_cols;
    return error;
}
}
