#include "dimensionreduction.h"
#include "dggsvd.h"
namespace DimReduction
{
/*
 * Structure Preserving Dimension Reduction for
 * Clustered Text Data Based on the Generalized Singular
 * Value Decomposition
 */
void lda_gsvd(
        const arma::mat& Hb,
        const arma::mat& Hw,
        arma::mat& G
        )
{
    arma::mat U,V,C,S,X;
    arma::mat A = Hb.t();
    arma::mat B = Hw.t();
    ML_Math::dggsvd(A,B,U,V,C,S,X);
    unsigned int N = G.n_cols;
    G = X.cols(0,N-1);
}
}
