#include "pcaplaneequ.h"
pcaplaneequ::pcaplaneequ(void)
{
    clear();
}

int pcaplaneequ::push_point(arma::fvec& p)// if the size of point is bigger than 3, we could calculate the equation, otherwise not.
{
    cloud.push_back(p);
    n++;

    ave = ave*float(n-1)/float(n)+cloud[n-1]/float(n);

    for (int i=0; i<n; i++)
    {
        cov_(0,0) = cov_(0,0) + (cloud[i](0)-ave(0))*(cloud[i](0)-ave(0))/float(n);
        cov_(0,1) = cov_(0,1) + (cloud[i](0)-ave(0))*(cloud[i](1)-ave(1))/float(n);
        cov_(0,2) = cov_(0,2) + (cloud[i](0)-ave(0))*(cloud[i](2)-ave(2))/float(n);
        cov_(1,1) = cov_(1,1) + (cloud[i](1)-ave(1))*(cloud[i](1)-ave(1))/float(n);
        cov_(1,2) = cov_(1,2) + (cloud[i](1)-ave(1))*(cloud[i](2)-ave(2))/float(n);
        cov_(2,2) = cov_(2,2) + (cloud[i](2)-ave(2))*(cloud[i](2)-ave(2))/float(n);
    }
    cov_(1,0) = cov_(0,1);
    cov_(2,0) = cov_(0,2);
    cov_(2,1) = cov_(1,2);

    arma::eig_sym(eval_,evec_,cov_,"std");

    int k=0;
    for (int i=1; i<3; i++)
    if (eval_(i)<eval_(k))
        k = i;
    a = evec_(0,k);
    b = evec_(1,k);
    c = evec_(2,k);
    d = -arma::dot(evec_.col(k),ave);
    eigennormal = eval_(k);
    normal = evec_.col(k);
    if (n<3)
        return 0;
    else
        return 1;
}

pcaplaneequ::~pcaplaneequ(void)
{
}
