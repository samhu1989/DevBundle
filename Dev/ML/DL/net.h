#ifndef NET_H
#define NET_H
#include "dl_global.h"
#include "optimizationcore.h"
#include "common.h"
#include <memory>
#include <layer.h>
namespace DL{
class DLSHARED_EXPORT Net:public Optimization::EnergyFunction
{
public:
    Net();
    virtual bool config(Config::Ptr);
    virtual inline size_t size(){return n_;}
    virtual void initialValue(arma::vec& x);
    virtual double gradient(const arma::vec& x,arma::vec& dx);
    virtual ~Net(){}
protected:
    virtual void forward();
    virtual void backward();
protected:
    std::vector<Layer::Ptr> layers_;
    size_t n_;
    double fx_;
    double* x_ptr_;
    double* dx_ptr_;
};
}
#endif // NET_H
