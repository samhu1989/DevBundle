#ifndef DIMENSIONREDUCTION_H
#define DIMENSIONREDUCTION_H
#include <armadillo>
#include "dimensionreduction_global.h"

namespace DimReduction
{
//comput lda using generalized singular value decomposition
void DIMENSIONREDUCTIONSHARED_EXPORT lda_gsvd(
        const arma::mat& Hb,
        const arma::mat& Hw,
        arma::mat& G
        );
}

#endif // DIMENSIONREDUCTION_H
