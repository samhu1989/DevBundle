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
    lapack_int m,n,p,l,k;
    m = A_.n_rows;
    n = A_.n_cols;
    p = B_.n_rows;
    k = 0;
    l = 0;
    arma::Col<lapack_int> iwork(n,arma::fill::zeros);
    arma::mat Q(n,n,arma::fill::zeros);
    arma::vec alpha(n,arma::fill::zeros);
    arma::vec beta(n,arma::fill::zeros);
    U = arma::mat(m,m,arma::fill::zeros);
    V = arma::mat(p,p,arma::fill::zeros);
    LAPACKE_dggsvd3(LAPACK_COL_MAJOR,
                   'U','V','Q',
                   m,n,p,&k,&l,
                   A_.memptr(),m,
                   B_.memptr(),p,
                   alpha.memptr(),
                   beta.memptr(),
                   U.memptr(),U.n_rows,
                   V.memptr(),V.n_rows,
                   Q.memptr(),Q.n_rows,
                   iwork.memptr()
                   );
    C = arma::mat(m,k+l,arma::fill::zeros);
    S = arma::mat(p,k+l,arma::fill::zeros);
    arma::mat R(k+l,k+l,arma::fill::zeros);
    arma::uword diag_len = 0;
    std::cerr<<"m:"<<m<<std::endl;
    std::cerr<<"n:"<<n<<std::endl;
    std::cerr<<"k:"<<k<<std::endl;
    std::cerr<<"l:"<<l<<std::endl;
    std::cerr<<iwork<<std::endl;
    if( m-k-l >= 0 )
    {
        diag_len = k+l;
//        std::cerr<<"1"<<std::endl;
        R = A_.submat(0,n-k-l,k+l-1,n-1);
    }else{
        diag_len = m;
//        std::cerr<<"2"<<std::endl;
        R.submat(0,0,m-1,k+l-1) = A_.submat(0,n-k-l,m-1,n-1);
        R.submat(m,m,k+l-1,k+l-1) = B_.submat(m-k,n+m-k-l,l-1,n-1);
    }
    C.diag().ones();
    arma::vec cd = C.diag();
    size_t Ncd = std::min( diag_len - k , cd.size() );
    cd.tail(Ncd) = alpha.subvec( k , Ncd + k - 1);
    C.diag() = cd;

    S.diag(1).ones();
    arma::vec sd = S.diag(1);
    size_t Nsd = std::min( diag_len - k , sd.size() );
    sd.head(Nsd) = beta.subvec( k , Nsd + k - 1);
    S.diag(1) = sd;

//    std::cerr<<"3"<<std::endl;
    arma::mat T(Q.n_rows,Q.n_cols,arma::fill::eye);
    T.submat(Q.n_rows-R.n_rows,Q.n_cols-R.n_cols,Q.n_rows-1,Q.n_cols-1) = R.i();
    X = Q*T;
}
}
#endif
