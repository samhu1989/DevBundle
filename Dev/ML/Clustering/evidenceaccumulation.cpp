#include "evidenceaccumulation.h"
namespace Clustering{
PEAC::PEAC()
{
    ;
}

bool PEAC::configure(Config::Ptr)
{
    return false;
}

void PEAC::compute(
        const arma::umat& E,
        const arma::umat& P,
        arma::uvec& y
        )
{
    iE_ = E;
    iP_ = P;
    compute();
    projectY(y);
}

void PEAC::compute()
{
    initY();
    computeCandN();
    initPList();
    initA();
    initG();
    initPQ();
    arma::uword t = 0 ;
    while( t < max_iter_ )
    {
        bestDY_ = pq_.front().value_;
        computeStep();
        computeY();
        computeA();
        computeD();
        computeG();
        updatePrior();
        ++t;
    }
}

void PEAC::computeCandN()
{
     N_ = iE_.n_rows;
     arma::uword index;
     #pragma omp parallel for
     for(index=0;index<iP_.n_cols;++index)
     {
         arma::uvec pair = iP_.col(index);
         arma::uvec sig0 = iE_.col(pair(0));
         arma::uvec sig1 = iE_.col(pair(1));
         arma::uvec index = arma::find( sig0 == sig1 );
         C_(pair(0),pair(1)) = double(index.size())/N_;
         C_(pair(1),pair(0)) = C_(pair(0),pair(1));
     }
}

void PEAC::initPList()
{
    ;
}

void PEAC::initY()
{
    ;
}

void PEAC::initA()
{
    ;
}

void PEAC::initG()
{
    ;
}

void PEAC::initPQ()
{
    ;
}

void PEAC::computeD()
{
    ;
}

void PEAC::computeStep()
{
    ;
}

void PEAC::computeY()
{
    ;
}

void PEAC::computeA()
{
    ;
}

void PEAC::computeG()
{
    ;
}

void PEAC::updatePrior()
{
    ;
}

void PEAC::projectY(arma::uvec& y)
{
    ;
}
}
