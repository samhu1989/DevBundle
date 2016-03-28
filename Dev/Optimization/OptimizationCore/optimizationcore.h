#ifndef OPTIMIZATIONCORE_H
#define OPTIMIZATIONCORE_H
#include <armadillo>
#include "optimizationcore_global.h"
namespace Optimization
{
class EnergyFunction{
public:
    virtual arma::vec initialValue() = 0;
    virtual double gradient( const arma::vec & x, arma::vec & dx ) = 0;
    static arma::vec numericGradient( EnergyFunction & efun, const arma::vec & x, double EPS=1e-3 );
    static arma::vec gradient( EnergyFunction & efun, const arma::vec & x );
    static double gradCheck( EnergyFunction & efun, const arma::vec & x, double EPS=1e-3 );
    static arma::fvec computeFunction( EnergyFunction & efun, const arma::vec & x, const arma::vec & dx, int n_samples = 100 );
    virtual ~EnergyFunction(){}
};

class Optimizer{
public:
    arma::vec minimize(EnergyFunction& efun,int restart ,bool verbose = false);
};
}
#endif // OPTIMIZATIONCORE_H
