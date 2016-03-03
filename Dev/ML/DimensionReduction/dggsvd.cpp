#ifndef DGGSVD_CPP
#define DGGSVD_CPP
#include "dggsvd.h"
#include <armadillo>
#include "OpenBLAS/lapacke.h"
#include "assert.h"
namespace ML_Math
{
void dggsvd(
        const arma::mat& A,
        const arma::mat& B,
        arma::mat& U,
        arma::mat& V,
        arma::mat& C,
        arma::mat& S,
        arma::mat& X
        )
{
    assert(A.n_cols==B.n_cols);
    arma::mat A_ = A;
    arma::mat B_ = B;
    lapack_int m,n,p,l,k,iwork;
    m = A_.n_rows;
    n = A_.n_cols;
    p = B_.n_rows;
    arma::mat Q(n,n,arma::fill::zeros);
    arma::vec alpha(n,arma::fill::zeros);
    arma::vec beta(n,arma::fill::zeros);
    U = arma::mat(m,m,arma::fill::zeros);
    V = arma::mat(p,p,arma::fill::zeros);
    LAPACKE_dggsvd(LAPACK_COL_MAJOR,
                   'U','V','Q',
                   m,n,p,&k,&l,
                   A_.memptr(),A_.n_rows,
                   B_.memptr(),B_.n_rows,
                   alpha.memptr(),
                   beta.memptr(),
                   U.memptr(),U.n_rows,
                   V.memptr(),V.n_rows,
                   Q.memptr(),Q.n_rows,
                   &iwork
                   );
    C = arma::mat(m,k+l,arma::fill::zeros);
    S = arma::mat(p,k+l,arma::fill::zeros);
    arma::mat R(k+l,k+l,arma::fill::zeros);
    arma::uword diag_len = 0;
    if( m-k-l >= 0 )
    {
        diag_len = k+l;
        R = A_.submat(0,n-k-l,k+l-1,n-1);
    }else{
        diag_len = m;
        R.submat(0,0,m-1,k+l-1) = A_.submat(0,n-k-l,m-1,n-1);
        R.submat(m-k,m,l-1,k+l-1) = B_.submat(m-k,n+m-k-l,l-1,n-1);
    }
    C.diag().ones();
    arma::vec cd = C.diag();
    cd.tail(diag_len - k) = alpha.subvec(k,diag_len-1);
    C.diag() = cd;

    S.diag(1).ones();
    arma::vec sd = S.diag(1);
    sd.head(diag_len - k) = beta.subvec(k,diag_len-1);
    S.diag(1) = sd;

    arma::mat T(Q.n_rows,Q.n_cols,arma::fill::eye);
    std::cerr<<"R:"<<std::endl;
    std::cerr<<R<<std::endl;
    std::cerr<<"T:"<<std::endl;
    std::cerr<<T<<std::endl;
    T.submat(Q.n_rows-R.n_rows,Q.n_cols-R.n_cols,Q.n_rows-1,Q.n_cols-1) = R.i();
    X = Q*T;
}
}
#endif
