#include "RegistrationBase.h"
#include "common.h"
#include "nanoflann.hpp"
namespace Registration{
void RegistrationBase::alignAroundZ(const arma::fmat &x0,const arma::fmat &x1, arma::fmat &bestR)
{
    arma::fmat::fixed<3,3> R;
    arma::fmat y;
    float N = 36.0;
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
    std::cerr<<"minError:"<<minError<<std::endl;
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
        error += dists(0);
    }
    error /= x1.n_cols;
    return error;
}
}
