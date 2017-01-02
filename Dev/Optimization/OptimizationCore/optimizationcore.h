#ifndef OPTIMIZATIONCORE_H
#define OPTIMIZATIONCORE_H
#include <armadillo>
#include "optimizationcore_global.h"
namespace Optimization
{
class OPTIMIZATIONCORESHARED_EXPORT EnergyFunction{
public:
    virtual size_t size()=0;
    virtual void initialValue(arma::vec& x) = 0;
    virtual double gradient( const arma::vec & x, arma::vec & dx ) = 0;
    static arma::vec numericGradient( EnergyFunction & efun, const arma::vec & x, float EPS=1e-3 );
    static arma::vec gradient( EnergyFunction & efun, const arma::vec & x );
    static double gradCheck( EnergyFunction & efun, const arma::vec & x, float EPS=1e-3 );
    static arma::vec computeFunction( EnergyFunction & efun, const arma::vec & x, const arma::vec & dx, int n_samples = 100 );
    virtual ~EnergyFunction(){}
};

class OPTIMIZATIONCORESHARED_EXPORT Optimizer{
public:
    void minimize(EnergyFunction& efun,arma::vec& x,int restart ,bool verbose = false);
};
}
#endif // OPTIMIZATIONCORE_H
