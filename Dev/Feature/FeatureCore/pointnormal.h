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
}
#include "pointnormal.hpp"
#endif // POINTNORMAL_H
