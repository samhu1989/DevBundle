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
    std::copy(x.memptr(),x.memptr()+x.size(),x_ptr_);
    forward();
    backward();
    std::copy(dx_ptr_,dx_ptr_+n_,dx.memptr());
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
