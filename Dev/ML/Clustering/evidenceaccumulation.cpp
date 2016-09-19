#include "evidenceaccumulation.h"
#include <assert.h>
#include <QTime>
namespace Clustering{
PEAC::PEAC():distr_(0.0001,1.0),rand_engine_(QTime::currentTime().msec())
{
    tol_ = std::numeric_limits<double>::epsilon();
    max_iter_ = 1000;
    k_ = 5;
}

bool PEAC::configure(Config::Ptr config)
{
    if(config->has("PEAC_Max_Iter"))
    {
        max_iter_ = config->getInt("PEAC_Max_Iter");
    }
    if(config->has("PEAC_tol"))
    {
        tol_ = config->getDouble("PEAC_tol");
    }
    if(config->has("PEAC_k"))
    {
        k_ = config->getInt("PEAC_k");
    }
    if(config->has("PEAC_init_mode"))
    {
        init_=config->getInt("PEAC_init_mode");
    }
    return true;
}

void PEAC::compute(
        const arma::umat& E,
        const arma::umat& P,
        arma::uvec& y
        )
{
    iE_ = E;
    iP_ = P;
//    std::cerr<<"PEAC::compute()"<<std::endl;
    compute();
//    std::cerr<<"projectY(y)"<<std::endl;
    projectY(y);
}

void PEAC::compute()
{
    computeCandN();
    if(init_==0)initY_RandIndex();
    if(init_==1)initY_Random();
    if(init_==2)initY_BestIndex();
    else last_obj_ = std::numeric_limits<double>::max();
    initA();
    initG();
    initPQ();
    arma::uword t = 0 ;
    double obj;
    while( t < max_iter_ )
    {
        std::cerr<<"bestDY=("<<bestDY_.alpha_<<","<<bestDY_.beta_<<","<<bestDY_.gamma_<<") with G:"<<pq_.front().prior_<<std::endl;
        if(pq_.front().prior_<0)break;
        if(pq_.front().prior_<tol_)break;
        getBestD();
        computeStep();
        computeY();
        computeA();
        computeG();
        updatePrior();
        obj = computeObj();
        std::cerr<<"obj("<<t<<")="<<obj<<std::endl;
        if(obj>last_obj_){
            *newY_ = *oldY_;
            break;
        }
        else last_obj_ = obj;
        ++t;
    }
}

void PEAC::computeCandN()
{
     N_ = iE_.n_rows;
     C_ = arma::sp_mat(iE_.n_cols,iE_.n_cols);
     #pragma omp parallel for
     for(arma::uword index=0;index<iP_.n_cols;++index)
     {
         arma::uvec pair = iP_.col(index);
         arma::uvec sig0 = iE_.col(pair(0));
         arma::uvec sig1 = iE_.col(pair(1));
         arma::uvec equal = arma::find( sig0 == sig1 );
         C_(pair(0),pair(1)) = double(equal.size())/double(N_);
         C_(pair(1),pair(0)) = C_(pair(0),pair(1));
     }
}

void PEAC::initY_Random()
{
    oldY_.reset(new arma::mat(k_,iE_.n_cols));
    oldY_->imbue([&](){ return distr_(rand_engine_);});
    arma::rowvec sum = arma::sum(*oldY_);
    oldY_->each_row() /= sum;

    newY_.reset(new arma::mat(k_,iE_.n_cols));
    *newY_ = *oldY_;
}

void PEAC::initY_RandIndex()
{
    std::uniform_int_distribution<arma::uword> rand(0,iE_.n_rows-1);
    arma::urowvec label = iE_.row(rand(rand_engine_));
    arma::uword max = arma::max(label);
    arma::uword min = arma::min(label);
    arma::uvec th = arma::linspace<arma::uvec>(min,max,k_+1);
    oldY_.reset(new arma::mat(k_,iE_.n_cols));
    newY_.reset(new arma::mat(k_,iE_.n_cols,arma::fill::zeros));
    #pragma omp parallel for
    for(arma::uword index=0;index<label.size();++index)
    {
        arma::uword k;
        for( k = th.size() - 1 ; k > 0 ;--k)
        {
            if( label(index) > th(k) )break;
        }
        assert( k >= 0);
        assert( k < newY_->n_rows );
        (*newY_)(k,index) = 1.0;
    }
    (*oldY_) = (*newY_);
}

void PEAC::initY_BestIndex()
{
    oldY_.reset(new arma::mat(k_,iE_.n_cols));
    newY_.reset(new arma::mat(k_,iE_.n_cols,arma::fill::zeros));
    double best_obj = std::numeric_limits<double>::max();
    arma::uword best_index = 0;
    double obj;
    for(arma::uword index=0;index<iE_.n_rows;++index)
    {
        arma::urowvec label = iE_.row(index);
        arma::uword max = arma::max(label);
        arma::uword min = arma::min(label);
        arma::uvec th = arma::linspace<arma::uvec>(min,max,k_+1);
        newY_->fill(0);
        #pragma omp parallel for
        for(arma::uword index=0;index<label.size();++index)
        {
            arma::uword k;
            for( k = th.size() - 1 ; k > 0 ;--k)
            {
                if( label(index) > th(k) )break;
            }
            assert( k >= 0);
            assert( k < newY_->n_rows );
            (*newY_)(k,index) = 1.0;
        }
        obj = computeObj();
//        std::cerr<<"Index:"<<index<<" with obj:"<<obj<<std::endl;
        if(obj<best_obj)
        {
            best_obj = obj;
            best_index = index;
        }
    }
    std::cerr<<"Best Index:"<<best_index<<" with obj:"<<best_obj<<std::endl;
    arma::urowvec label = iE_.row(best_index);
    arma::uword max = arma::max(label);
    arma::uword min = arma::min(label);
    arma::uvec th = arma::linspace<arma::uvec>(min,max,k_+1);
    newY_->fill(0);
    #pragma omp parallel for
    for(arma::uword index=0;index<label.size();++index)
    {
        arma::uword k;
        for( k = th.size() - 1 ; k > 0 ;--k)
        {
            if( label(index) > th(k) )break;
        }
        assert( k >= 0);
        assert( k < newY_->n_rows );
        (*newY_)(k,index) = 1.0;
    }
    last_obj_ = best_obj;
    (*oldY_) = (*newY_);
}

void PEAC::initA()
{
    oldA_.reset(new arma::sp_mat(iE_.n_cols,iE_.n_cols));
    newA_.reset(new arma::sp_mat(iE_.n_cols,iE_.n_cols));
    for(arma::uword index=0;index<iP_.n_cols;++index)
    {
        arma::uvec pair = iP_.col(index);
        arma::uword i = pair(0);
        arma::uword j = pair(1);
        (*oldA_)(i,j) = arma::accu(newY_->col(i)%newY_->col(j));
        (*newA_)(i,j) = (*oldA_)(i,j);
    }
}

void PEAC::initG()
{
    G_ = arma::mat(k_,iE_.n_cols,arma::fill::zeros);
    for(arma::uword index=0;index<iP_.n_cols;++index)
    {
        arma::uvec pair = iP_.col(index);
        arma::uword i = pair(0);
        arma::uword j = pair(1);
        G_.col(i) += N_*newY_->col(j)*((*newA_)(i,j) - C_(i,j));
        G_.col(j) += N_*newY_->col(i)*((*newA_)(i,j) - C_(i,j));
    }
}

void PEAC::initPQ()
{
    pq_.resize(iE_.n_cols);
    #pragma omp parallel for
    for(arma::uword i=0;i<iE_.n_cols;++i)
    {
        Triplet& tri = pq_[i];
        tri.index_ = i;
        computePrior(tri);
    }
    std::make_heap(pq_.begin(),pq_.end(),std::less<Triplet>());
}

void PEAC::getBestD()
{
    bestDY_ = pq_.front().value_;
    std::vector<arma::uword> Pgamma_vec;
    Pgamma_vec.reserve(iE_.n_cols);
    for(arma::uword index=0;index<iP_.n_cols;++index)
    {
        arma::uvec pair = iP_.col(index);
        arma::uword i = pair(0);
        arma::uword j = pair(1);
        if(i==bestDY_.gamma_)
        {
            Pgamma_vec.push_back(j);
        }
        if(j==bestDY_.gamma_)
        {
            Pgamma_vec.push_back(i);
        }
    }
    Pgamma_ = arma::uvec(Pgamma_vec);
}

void PEAC::computeStep()
{
    double d;
    d = G_(bestDY_.beta_,bestDY_.gamma_) - G_(bestDY_.alpha_,bestDY_.gamma_);
    arma::rowvec Ya = (*newY_).row(bestDY_.alpha_);
    arma::rowvec Yb = (*newY_).row(bestDY_.beta_);
    double frac = arma::accu(arma::square(Ya(Pgamma_) - Yb(Pgamma_)));
    d /= ( N_*frac );
    step_ = std::min(d,(*newY_)(bestDY_.beta_,bestDY_.gamma_));
}

void PEAC::computeY()
{
    oldY_.swap(newY_);
    *newY_ = *oldY_;
    (*newY_)(bestDY_.alpha_,bestDY_.gamma_) += step_;
    (*newY_)(bestDY_.beta_,bestDY_.gamma_) -= step_;
}

void PEAC::computeA()
{
    oldA_.swap(newA_);
    *newA_ = *oldA_;
    arma::rowvec Ya = (*oldY_).row(bestDY_.alpha_);
    arma::rowvec Yb = (*oldY_).row(bestDY_.beta_);
    arma::rowvec dif = step_*(Ya - Yb);
    for(arma::uword index=0;index<Pgamma_.n_rows;++index)
    {
        arma::uword i = bestDY_.gamma_;
        arma::uword j = Pgamma_(index);
        (*newA_)(i,j) += dif(j);
        (*newA_)(i,j) = (*newA_)(j,i);
    }
}

void PEAC::computeG()
{
    //if j \in P_{\gamma}
    G_.col(bestDY_.gamma_).fill(0.0);
    for(arma::uvec::iterator iter=Pgamma_.begin();iter!=Pgamma_.end();++iter)
    {
        arma::uword j = *iter;
        arma::vec a = (*newY_).col(bestDY_.gamma_)*( (*newA_)(bestDY_.gamma_,j) - C_(bestDY_.gamma_,j) );
        arma::vec b = (*oldY_).col(bestDY_.gamma_)*( (*oldA_)(bestDY_.gamma_,j) - C_(bestDY_.gamma_,j) );
        G_.col(j) += N_*( a - b );
        G_.col(bestDY_.gamma_) += N_*(*newY_).col(j)*( (*newA_)(j,bestDY_.gamma_) - C_(j,bestDY_.gamma_) );
    }
}

void PEAC::computePrior(Triplet& tri)
{
    arma::uword i = tri.index_;
    arma::vec gY = G_.col(i);
    gY.min(tri.value_.alpha_);
    arma::uvec index = arma::find(newY_->col(i)!=0);
    arma::vec subgY = gY(index);
    arma::uword mi;
    subgY.max(mi);
    tri.value_.beta_ = index(mi);
    tri.value_.gamma_ = i;
    tri.prior_ = gY(tri.value_.beta_) - gY(tri.value_.alpha_);
}

double PEAC::computeObj()
{
    double obj = 0.0;
    for(arma::uword index=0;index<iP_.n_cols;++index)
    {
        arma::uvec pair = iP_.col(index);
        arma::uword i = pair(0);
        arma::uword j = pair(1);
        double dif = C_(i,j) - arma::accu((*newY_).col(i)%(*newY_).col(j));
        obj += N_*(dif*dif);
    }
    return obj;
}

void PEAC::updatePrior()
{
    arma::uvec indexed(pq_.size());
    #pragma omp parallel for
    for(arma::uword i=0;i<pq_.size();++i)
    {
        indexed(pq_[i].index_)=i;
    }
    computePrior(pq_[indexed[bestDY_.gamma_]]);
    #pragma omp parallel for
    for(arma::uword i=0;i<Pgamma_.size();++i)
    {
        computePrior(pq_[indexed[Pgamma_[i]]]);
    }
    std::make_heap(pq_.begin(),pq_.end(),std::less<Triplet>());
}

void PEAC::projectY(arma::uvec& y)
{
    y = arma::uvec(iE_.n_cols);
    #pragma omp parallel for
    for(arma::uword index=0;index<newY_->n_cols;++index)
    {
        arma::uword mi;
        newY_->col(index).max(mi);
        y(index) = mi + 1;
    }
}
}
