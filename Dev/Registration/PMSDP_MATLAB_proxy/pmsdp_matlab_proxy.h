#ifndef PMSDP_MATLAB_PROXY_H
#define PMSDP_MATLAB_PROXY_H

#include "pmsdp_matlab_proxy_global.h"
#include <armadillo>
extern "C"{
bool PMSDP_MATLAB_PROXYSHARED_EXPORT initMatlab(void);
void PMSDP_MATLAB_PROXYSHARED_EXPORT compute(
        const arma::mat &P,
        const arma::mat &Q,
        arma::mat& R,
        arma::uvec& X
);
void PMSDP_MATLAB_PROXYSHARED_EXPORT terminateMatlab(void);
}
#endif // PMSDP_MATLAB_PROXY_H
