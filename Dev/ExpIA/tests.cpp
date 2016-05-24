#ifndef TESTS_CPP
#define TESTS_CPP
#include "tests.h"
#include "OpenBLAS/lapacke.h"
#include <armadillo>
#include <assert.h>
#include "dggsvd.h"
#include "common.h"
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
void Inside_BBox_test()
{
    bool sucess = true;
    arma::fmat box={{0,1,1,0},{0,0,1,1}};
    arma::fvec p1={0.5,0.5};
    if(!in2DMBB(box,p1))
    {
        sucess = false;
        std::cerr<<"Failed for"<<std::endl;
        std::cerr<<box<<std::endl;
        std::cerr<<p1<<std::endl;
    }
    arma::fvec p2={0.5,1.5};
    if(in2DMBB(box,p2))
    {
        sucess = false;
        std::cerr<<"Failed for"<<std::endl;
        std::cerr<<box<<std::endl;
        std::cerr<<p2<<std::endl;
    }
    arma::fvec p3={-0.5,-0.5};
    if(in2DMBB(box,p3))
    {
        sucess = false;
        std::cerr<<"Failed for"<<std::endl;
        std::cerr<<box<<std::endl;
        std::cerr<<p3<<std::endl;
    }
    arma::fvec p4={0.3,0.4};
    if(!in2DMBB(box,p4))
    {
        sucess = false;
        std::cerr<<"Failed for"<<std::endl;
        std::cerr<<box<<std::endl;
        std::cerr<<p4<<std::endl;
    }
    if(sucess)
    {
        std::cerr<<"Success in Inside_BBox_test()"<<std::endl;
    }
}
}
#endif // TESTS
