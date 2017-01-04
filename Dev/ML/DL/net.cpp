#include "net.h"
#include <cassert>
namespace DL{
Net::Net()
{
    ;
}

bool Net::config(Config::Ptr)
{
    ;
}

void Net::initialValue(arma::vec& x)
{
    if((double*)x.memptr()!=x_ptr_ || x.size() != n_ )
    {
        x = arma::vec(x_ptr_,n_,false,true);
    }
}

double Net::gradient(const arma::vec& x,arma::vec& dx)
{
    assert(x.size()!=n_);
    assert(dx.size()!=n_);
    x_ptr_ = (double*)x.memptr();
    dx_ptr_ = dx.memptr();
    forward();
    backward();
    return fx_; // return function value
}

void Net::forward()
{
    ;
}

void Net::backward()
{
    ;
}

}
