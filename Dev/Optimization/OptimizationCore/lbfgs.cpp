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
    arma::rowvec tmp((double*)x,n,false,true);
    std::cerr<<"x:"<<tmp<<std::endl;
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
    arma::vec vg((lbfgsfloatval_t*)g,n,false,true);
    arma::vec vx((lbfgsfloatval_t*)x,n,false,true);
//    arma::fvec vx( n ), vg( n );
//    std::copy( x, x+n, (lbfgsfloatval_t*)( vx.memptr()) );
    lbfgsfloatval_t r = efun->gradient( vx, vg );
//    std::copy( (lbfgsfloatval_t*)vg.memptr(), (lbfgsfloatval_t*)(vg.memptr())+n, g );
    return r;
}

void LBFGS::minimize(EnergyFunction& efun, arma::vec &x0, int restart, bool verbose)
{
    const int n = efun.size();
    lbfgsfloatval_t *x = lbfgs_malloc(n);

    if( x0.size() != n || !x0.is_finite() ){
        x0 = arma::vec(x,n,false,true);
        efun.initialValue(x0);
    }else {
        std::copy((lbfgsfloatval_t*)x0.memptr(),((lbfgsfloatval_t*)x0.memptr())+n,x);
        x0 = arma::vec(x,n,false,true);
        efun.initialValue(x0);
    }

    if (x == NULL) {
        printf("ERROR: Failed to allocate a memory block for variables.\n");
        return;
    }

    lbfgs_parameter_t param;
    lbfgs_parameter_init(&param);
    // You might want to adjust the parameters to your problem
    param.epsilon = 1e-6;
    param.max_iterations = 50;

    float last_f = 1e100;
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
        if(0!=ret)
        {
            printf("LBFGS_SUCCESS:%d\n",LBFGS_SUCCESS);
            printf("LBFGS_CONVERGENCE:%d\n",LBFGS_CONVERGENCE);
            printf("LBFGS_STOP:%d\n",LBFGS_STOP);
            printf("LBFGS_ALREADY_MINIMIZED:%d\n",LBFGS_ALREADY_MINIMIZED);
            printf("LBFGSERR_UNKNOWNERROR:%d\n",LBFGSERR_UNKNOWNERROR);
            printf("LBFGSERR_LOGICERROR:%d\n",LBFGSERR_LOGICERROR);
            printf("LBFGSERR_OUTOFMEMORY:%d\n",LBFGSERR_OUTOFMEMORY);
            printf("LBFGSERR_CANCELED:%d\n",LBFGSERR_CANCELED);
            printf("LBFGSERR_INVALID_N:%d\n",LBFGSERR_INVALID_N);
            printf("LBFGSERR_INVALID_N_SSE:%d\n",LBFGSERR_INVALID_N_SSE);
            printf("LBFGSERR_INVALID_X_SSE:%d\n",LBFGSERR_INVALID_X_SSE);
            printf("LBFGSERR_INVALID_EPSILON:%d\n",LBFGSERR_INVALID_EPSILON);
            printf("LBFGSERR_INVALID_TESTPERIOD:%d\n",LBFGSERR_INVALID_TESTPERIOD);
            printf("LBFGSERR_INVALID_DELTA:%d\n",LBFGSERR_INVALID_DELTA);
            printf("LBFGSERR_INVALID_LINESEARCH:%d\n",LBFGSERR_INVALID_LINESEARCH);
            printf("LBFGSERR_INVALID_MINSTEP:%d\n",LBFGSERR_INVALID_MINSTEP);
            printf("LBFGSERR_INVALID_MAXSTEP:%d\n",LBFGSERR_INVALID_MAXSTEP);
            printf("LBFGSERR_INVALID_FTOL:%d\n",LBFGSERR_INVALID_FTOL);
            printf("LBFGSERR_INVALID_WOLFE:%d\n",LBFGSERR_INVALID_WOLFE);
            printf("LBFGSERR_INVALID_GTOL:%d\n",LBFGSERR_INVALID_GTOL);
            printf("LBFGSERR_INVALID_XTOL:%d\n",LBFGSERR_INVALID_XTOL);
            printf("LBFGSERR_INVALID_MAXLINESEARCH:%d\n",LBFGSERR_INVALID_MAXLINESEARCH);
            printf("LBFGSERR_INVALID_ORTHANTWISE:%d\n",LBFGSERR_INVALID_ORTHANTWISE);
            printf("LBFGSERR_INVALID_ORTHANTWISE_START:%d\n",LBFGSERR_INVALID_ORTHANTWISE_START);
            printf("LBFGSERR_INVALID_ORTHANTWISE_END:%d\n",LBFGSERR_INVALID_ORTHANTWISE_END);
            printf("LBFGSERR_OUTOFINTERVAL:%d\n",LBFGSERR_OUTOFINTERVAL);
            printf("LBFGSERR_INCORRECT_TMINMAX:%d\n",LBFGSERR_INCORRECT_TMINMAX);
            printf("LBFGSERR_ROUNDING_ERROR:%d\n",LBFGSERR_ROUNDING_ERROR);
            printf("LBFGSERR_MINIMUMSTEP:%d\n",LBFGSERR_MINIMUMSTEP);
            printf("LBFGSERR_MAXIMUMSTEP:%d\n",LBFGSERR_MAXIMUMSTEP);
            printf("LBFGSERR_MAXIMUMLINESEARCH:%d\n",LBFGSERR_MAXIMUMLINESEARCH);
            printf("LBFGSERR_MAXIMUMITERATION:%d\n",LBFGSERR_MAXIMUMITERATION);
            printf("LBFGSERR_INVALIDPARAMETERS:%d\n",LBFGSERR_WIDTHTOOSMALL);
            printf("LBFGSERR_INVALIDPARAMETERS:%d\n",LBFGSERR_INVALIDPARAMETERS);
            printf("LBFGSERR_INCREASEGRADIENT:%d\n",LBFGSERR_INCREASEGRADIENT);
        }
    }

    std::copy( x, x+n, (double*)(x0.memptr()) );
    lbfgs_free(x);
    return;
}
}
