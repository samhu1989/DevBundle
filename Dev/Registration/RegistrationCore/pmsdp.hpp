#ifndef PMSDP_HPP
#define PMSDP_HPP
#include "pmsdp.h"

#include <QStandardPaths>
#include <QDebug>
namespace Registration {
template<typename M>
PMSDP<M>::PMSDP()
{

}
template<typename M>
PMSDP<M>::~PMSDP()
{

}

template<typename M>
bool PMSDP<M>::configure(Config::Ptr& config,InfoPtr& info)
{
    info = std::make_shared<Info>();
    if(!config)
    {
        return true;
    }
    return false;
}

template<typename M>
void PMSDP<M>::computeByMatlab(void)
{
    std::cerr<<"PMSDP::computeByMatlab(void)"<<std::endl;
    if(!pmsdp_matlab_){
        pmsdp_matlab_=std::make_shared<PMSDP_MATLAB>("./bin/PMSDP_MATLAB_proxy.dll");
        if(!pmsdp_matlab_->init())return;
    }
    arma::mat R;
    arma::uvec X;
    arma::mat P = arma::conv_to<arma::mat>::from(*P_);
    arma::mat Q = arma::conv_to<arma::mat>::from(*Q_);
    pmsdp_matlab_->compute(P,Q,R,X);
    *P_ = arma::conv_to<arma::fmat>::from(R*P);
}

template<typename M>
void PMSDP<M>::generateObj()
{
    std::cerr<<"PMSDP::generateObj()"<<std::endl;
    std::cerr<<"Under Implementation"<<std::endl;
//    arma::uword n = P_.n_cols;
//    arma::uword d = P_.n_rows;
//    arma::uword size = n+d*d+1;
//    std::vector<arma::mat>C(n,arma::mat(size,size,arma::fill::zeros));
//    arma::mat w = arma::kron(P_,Q_);
//    for(arma::uword ic=0;ic<n;++ic)
//    {
//        arma::mat& block = C[ic];
//        arma::mat wi = w.cols(ic*n,(ic+1)*n-1);
//        block.submat(1,n+1,n,n+d*d) = -1.0*wi.t();
//        block.submat(n+1,1,n+d*d,n) = -1.0*wi;
//    }
//    sdp_.setC(C);
}

template<typename M>
void PMSDP<M>::generateConstraint()
{
    std::cerr<<"PMSDP::generateConstraint()"<<std::endl;
    std::cerr<<"Under Implementation"<<std::endl;
//    std::vector<std::vector<arma::mat>> As;
//    std::vector<double> b;
    /*
     * A{ii} = diag(X(:,ii));
     */
//    sdp_.setAs(As);
//    sdp_.setb(b);
}
template<typename M>
void PMSDP<M>::projectRX()
{
    std::cerr<<"PMSDP::projectRX()"<<std::endl;
    std::cerr<<"Under Implementation"<<std::endl;
}

template<typename M>
bool PMSDP<M>::initForThread(void* meshlistptr,InfoPtr&info){
    MeshList* list = reinterpret_cast<MeshList*>(meshlistptr);
    if(!list){
        error_string_ = "Can not locate the inputs";
        return false;
    }
    if(list->size()!=2){
        error_string_ = "This algorithm is designed for two input meshes";
        return false;
    }
    MeshBundle<DefaultMesh>::Ptr m0 = (*list)[0];
    P_=std::make_shared<arma::fmat>(
             (float*)m0->mesh_.points(),
              3,
             m0->mesh_.n_vertices(),
             false,
             true
             );
    MeshBundle<DefaultMesh>::Ptr m1 = (*list)[1];
    Q_=std::make_shared<arma::fmat>(
                (float*)m1->mesh_.points(),
                3,
                m1->mesh_.n_vertices(),
                false,
                true
                );
    return true;
}
}
#endif // PMSDP_HPP
