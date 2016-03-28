#include "lbfgs.h"
#include <exception>
namespace Optimization{
LBFGS::LBFGS()
{

}
int LBFGS::progress(
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
        )
{
    std::cerr<<"Iteration "<<k<<std::endl;
    std::cerr<<"fx = "<<fx<<", xnorm = "<<xnorm<<", gnorm = "<<gnorm<<", step = "<<step<<std::endl;
    std::cerr<<std::endl;
    return 0;
}

lbfgsfloatval_t LBFGS::evaluate(
        void *instance,
        const lbfgsfloatval_t *x,
        lbfgsfloatval_t *g,
        const int n,
        const lbfgsfloatval_t step
        )
{

    EnergyFunction * efun = static_cast<EnergyFunction*>( instance );
    if(!efun)throw std::logic_error("invalid energy function instance");
    arma::vec vx( n ), vg( n );
    std::copy( x, x+n, (lbfgsfloatval_t*)( vx.memptr()) );
    lbfgsfloatval_t r = efun->gradient( vx, vg );
    std::copy( (lbfgsfloatval_t*)vg.memptr(), (lbfgsfloatval_t*)(vg.memptr())+n, g );
    return r;
}

arma::vec LBFGS::minimize(EnergyFunction& efun,int restart,bool verbose)
{
    arma::vec x0 = efun.initialValue();
    const int n = x0.n_rows;

    lbfgsfloatval_t *x = lbfgs_malloc(n);
    if (x == NULL) {
        printf("ERROR: Failed to allocate a memory block for variables.\n");
        return x0;
    }
    std::copy((lbfgsfloatval_t*)x0.memptr(),((lbfgsfloatval_t*)x0.memptr())+n,x);

    lbfgs_parameter_t param;
    lbfgs_parameter_init(&param);
    // You might want to adjust the parameters to your problem
    param.epsilon = 1e-6;
    param.max_iterations = 50;

    double last_f = 1e100;
    int ret;
    for( int i=0; i<=restart; i++ ) {
        lbfgsfloatval_t fx;
        ret = lbfgs(n, x, &fx, LBFGS::evaluate, verbose?LBFGS::progress:NULL, &efun, &param);
        if( last_f > fx )
            last_f = fx;
        else
            break;
    }

    if ( verbose ) {
        printf("L-BFGS optimization terminated with status code = %d\n", ret);
    }

    std::copy( x, x+n, (double*)(x0.memptr()) );
    lbfgs_free(x);
    return x0;
}
}
