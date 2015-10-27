#include "coherentpointdrift.h"
#include <cmath>
namespace Registration
{
    void CPDBase::stepE()
    {
        //initialize the probability matrix
        if(P_.is_empty())
        {
            P_ = arma::fmat(Y_.n_cols,X_.n_cols);
        }
        //by default D = 3
        float c = std::pow(2*M_PI*var,1.5)*(InfoPtr_->omega/(1-InfoPtr_->omega))*(float(Y_.n_cols)/float(X_.n_cols));
        float k = - 0.5 / var;
        for(int r=0;r<P_.n_rows;++r)
        {
            arma::fmat tmpm = X_.each_col() - Y_.col(r);
            P_.row(r) = arma::sum(arma::square(tmpm));
        }
        P_ = arma::exp(k*P_);
        arma::fvec p_sum = arma::sum(P_) + c;
        for(int c=0;c<P_.n_cols;++c)
        {
            P_.col(c) /= p_sum(c);
        }
    }


}


