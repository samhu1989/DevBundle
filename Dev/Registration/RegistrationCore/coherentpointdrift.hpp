#ifndef COHERENTPOINTDRIFT_HPP
#define COHERENTPOINTDRIFT_HPP
#include "coherentpointdrift.h"
namespace Registration{
template<typename M>
CPDRigid3D<M>::CPDRigid3D():CPDBase()
{
    ResPtr_ = std::shared_ptr<Result>(new Result);
}

template<typename M>
CPDRigid3D<M>::~CPDRigid3D()
{
    ;
}

template<typename M>
void CPDRigid3D<M>::stepM()
{
    //update Center
    dividebyNp_ = 1.0 / arma::accu(P_);
    muX = dividebyNp_*arma::sum(X_*(P_.t()),2);
    mX_ = X_.each_col() - muX;
    muY = dividebyNp_*arma::sum(Y_*P_,2);
    mY_ = Y_.each_col() - muY;

    //update Rotation
    A_ = mX_*P_.t()*mY_.t();
    arma::fmat U;
    arma::fmat V;
    arma::fvec s;
    arma::fmat R(ResPtr_->R,3,3,false,true);
    if(!arma::svd(U,s,V,A_,"std"))
    {
        std::cerr<<"Faild to do svd on A!"<<std::endl;
    }
    arma::fmat C(3,3,arma::fill::eye);
    C(2,2) = arma::det( U * V.t() )>=0 ? 1.0 : -1.0;
    R = U*C*(V.t());

    //trace( A_.t() * R ) with A_ and R are 3x3
    float trAtR = arma::accu(A_%R);

    //trace( mY_ * d(P1) * mY_ )
    arma::fvec P_SumRow_ = arma::sum(P_,2);
    arma::fmat square_mY_ = arma::square(mY_);
    square_mY_.each_row() %= P_SumRow_;
    float trYPY = arma::accu(square_mY_);

    //trace( mX.t() * d(P.t()1) * mX )
    arma::fvec P_SumCol_ = arma::sum(P_,1);
    arma::fmat square_mX_ = arma::square(mX_);
    square_mX_.each_row() %= P_SumCol_;
    float trXPX = arma::accu(square_mX_);

    //update Scale if required
    if(InfoPtr_->isScaled)
    {
        if( trAtR != 0.0 )ResPtr_->s = trAtR / trYPY;
        else{
            std::cerr<<"Get wrong b in updateScale"<<std::endl;
        }
    }

    //update Translation
    arma::fvec t(ResPtr_->t,3,false,true);
    t = muX - ResPtr_->s*(R*muY);

    //update Var
    var = dividebyNp_*( trXPX - ResPtr_->s*trAtR ) / 3.0 ;

    //update Y
    Y_ = R*Y_ ;
    Y_.each_col() += t;
}

template<typename M>
void CPDRigid3D<M>::reset(const arma::fmat&source,const arma::fmat&target,Info&info)
{
    CPDBase::reset(source,target,info);
    arma::fmat R(ResPtr_->R,3,3,false,true);
    R = arma::eye<arma::fmat>(3,3);
    arma::fvec t(ResPtr_->t,3,false,true);
    t = arma::zeros<arma::fmat>(3);
    ResPtr_->s = 1.0;
    info.result = (void*)ResPtr_.get();
}

template<typename M>
float CPDRigid3D<M>::fitness()
{
    ;
}

template<typename M>
bool CPDRigid3D<M>::isEnd()
{
    if( count >= InfoPtr_->max_iter ){
        InfoPtr_->mode = CPDBase::MaxIter;
        return true;
    }
    if( var < InfoPtr_->var_th ){
        InfoPtr_->mode = CPDBase::VarBelowThreshold;
        return true;
    }
    if( fitness() < InfoPtr_->fitness_th ){
        InfoPtr_->mode = CPDBase::ErrorBelowThreshold;
        return true;
    }
    return false;
}
template<typename M>
void CPDRigid3D<M>::configure(Info&)
{
    ;
}

template<typename M>
bool CPDRigid3D<M>::initForThread(void* meshlistptr)
{
    MeshList* list = reinterpret_cast<MeshList*>(meshlistptr);
    if(!list)return false;
    if(list->size()!=2){
        std::clog<<"This algorithm is designed for two input meshes"<<std::endl;
        return false;
    }
    arma::fmat source((float*)(*list)[0]->mesh_.points(),3,(*list)[0]->mesh_.n_vertices(),false,true);
    arma::fmat target((float*)(*list)[1]->mesh_.points(),3,(*list)[1]->mesh_.n_vertices(),false,true);

    std::shared_ptr<Info> info(new Info);
    configure(*info);

    reset(source,target,*info);
}

}
#endif // COHERENTPOINTDRIFT_HPP

