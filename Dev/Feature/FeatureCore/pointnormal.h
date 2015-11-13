#ifndef POINTNORMAL_H
#define POINTNORMAL_H
#include "common.h"
#include <armadillo>
#include "featurecore_global.h"
namespace Feature
{
template<typename M>
void computePointNormal(M& mesh,float r=0.1,int k=10);

template<typename M>
void computePointNormal(M& mesh,std::shared_ptr<float>& curvature,float r=0.1,int k=10);

inline void fitPlane(
        arma::fvec &center,
        arma::fmat &neighbor,
        arma::fvec &normal,
        float& distToOrigin)
{
    neighbor.each_col() -= center;
    arma::fmat c = arma::cov( neighbor.t() );
    arma::fmat U,V;
    arma::fvec s;
    if(!arma::svd(U,s,V,c,"std"))
    {
        normal.fill(std::numeric_limits<float>::quiet_NaN());
        distToOrigin = std::numeric_limits<float>::quiet_NaN();
    }else{
        arma::uword minidx;
        s.min(minidx);
        normal = U.col(minidx);
        distToOrigin = - arma::dot(normal,center);
    }
}

inline void fitPlane(
        arma::fvec &center,
        arma::fmat &neighbor,
        arma::fvec &normal,
        float &curvature,
        float& distToOrigin)
{
    neighbor.each_col() -= center;
    arma::fmat c = arma::cov( neighbor.t() );
    arma::fmat U,V;
    arma::fvec s;
    if(!arma::svd(U,s,V,c,"std"))
    {
        normal.fill(std::numeric_limits<float>::quiet_NaN());
        distToOrigin = std::numeric_limits<float>::quiet_NaN();
        curvature = std::numeric_limits<float>::quiet_NaN();
    }else{
        arma::uword minidx;
        curvature = s.min(minidx);
        normal = U.col(minidx);
        distToOrigin = - arma::dot(normal,center);
    }
}
}
#include "pointnormal.hpp"
#endif // POINTNORMAL_H
