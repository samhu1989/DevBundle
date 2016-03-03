#ifndef DGGSVD_H
#define DGGSVD_H
#include "dimensionreduction_global.h"
#include <armadillo>
namespace ML_Math
{
/*
 * U'AX = [0 C]
 * V'BX = [0 S]
 */
void DIMENSIONREDUCTIONSHARED_EXPORT dggsvd(
        const arma::mat& A,
        const arma::mat& B,
        arma::mat& U,
        arma::mat& V,
        arma::mat& C,
        arma::mat& S,
        arma::mat& X
        );
}
#endif // DGGSVD

