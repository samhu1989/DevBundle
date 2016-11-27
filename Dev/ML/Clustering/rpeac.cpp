#include "rpeac.h" 
namespace Clustering {
RPEAC::RPEAC():PEAC()
{

}
void RPEAC::convert_to_atomic(
        const arma::umat&E,
        const arma::umat&P,
        arma::umat &aE,
        arma::umat &aP
        )
{
    std::cerr<<"E("<<E.n_rows<<","<<E.n_cols<<")"<<std::endl;
    sig_hash_.clear();
    index_map_.resize(E.n_cols);
    for( arma::uword i=0 ; i < E.n_cols ; ++i )
    {
        //first appear
        if( sig_hash_.find( E.col(i) ) == sig_hash_.end() )
        {
            sig_hash_[ E.col(i) ] = sig_hash_.size() - 1;
        }
        index_map_[i] = sig_hash_[ E.col(i) ] ;
//        index_map_[i] = E(12,i);
    }
    aE = arma::umat(E.n_rows,sig_hash_.size());
    for (
         uvec_hash_map::iterator iter = sig_hash_.begin()   ;
                                 iter != sig_hash_.end()    ;
                                            ++iter
         )
    {

        aE.col(iter->second) = iter->first;
    }
    std::vector<arma::uword> P_vec;
    uvec_hash_map edge_hash_;
    arma::uvec pair(2);
    for(arma::umat::const_iterator iter=P.begin();iter!=P.end(); )
    {
        pair(0) = *iter;
        ++iter;
        pair(1) = *iter;
        ++iter;
        if( edge_hash_.find(pair) == edge_hash_.end() ) //first appear
        {
            edge_hash_[ pair ] = edge_hash_.size() - 1;
            P_vec.push_back(index_map_[pair(0)]);
            P_vec.push_back(index_map_[pair(1)]);
        }
    }
    aP = arma::umat(P_vec.data(),2,P_vec.size()/2,true,true);
}
void RPEAC::compute(
        const arma::umat& E,
        const arma::umat& P,
        arma::uvec& y
        )
{
    arma::umat aE;
    arma::umat aP;
    arma::uvec ay;
    convert_to_atomic(E,P,aE,aP);
    PEAC::compute(aE,aP,ay);
    std::cerr<<"atomic number:"<<ay.n_rows<<std::endl;
    y = arma::uvec(E.n_cols,arma::fill::zeros);
    convert_from_atomic(ay,y);
//    show_atomic(y);
}
void RPEAC::show_atomic(arma::uvec& y)
{
    y = arma::uvec(index_map_);
}

void RPEAC::convert_from_atomic(
        const arma::uvec& ay,
        arma::uvec& y
        )
{
    #pragma omp parallel for
    for( arma::uword i=0 ; i < y.n_rows ; ++i )
    {
        y(i) = ay(index_map_[i]);
    }
}
}
