#ifndef TESTS_CPP
#define TESTS_CPP
#include "tests.h"
#include "OpenBLAS/lapacke.h"
#include <armadillo>
#include <assert.h>
#include "dggsvd.h"
namespace TEST {
void LAPACKE_dggsvd_test()
{
    arma::mat A={{1,2,3},{3,2,1},{4,5,6},{7,8,8}};
    arma::mat B={{-2,-3,3},{4,6,5}};
    arma::mat U,V,C,S,X;
    ML_Math::dggsvd(A,B,U,V,C,S,X);
    std::cerr<<U*C*X.i()<<std::endl;
    std::cerr<<V*S*X.i()<<std::endl;
}
}
#endif // TESTS
