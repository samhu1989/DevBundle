#ifndef LBFGS_H
#define LBFGS_H
#include "optimizationcore.h"
#include "LBFGS/lbfgscore.h"
namespace Optimization{
class LBFGS : public Optimizer
{
public:
    LBFGS();
    virtual ~LBFGS(){}
    virtual arma::vec minimize(EnergyFunction& efun,int restart = 0,bool verbose = false);
    static lbfgsfloatval_t evaluate(
            void *instance,
            const lbfgsfloatval_t *x,
            lbfgsfloatval_t *g,
            const int n,
            const lbfgsfloatval_t step
            );
    static int progress(
            void *instance,
            const lbfgsfloatval_t *x,
            const lbfgsfloatval_t *g,
            const lbfgsfloatval_t fx,
            const lbfgsfloatval_t xnorm,
            const lbfgsfloatval_t gnorm,
            const lbfgsfloatval_t step,
            int n,
            int k,
            int ls
            );
};
}
#endif // LBFGS_H
