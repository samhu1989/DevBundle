#ifndef RPEAC_H
#define RPEAC_H
#include "evidenceaccumulation.h"
#include <ext/hash_map>
#include <hash_fun.h>
namespace Clustering{
class CLUSTERINGSHARED_EXPORT RPEAC : public PEAC
{
public:
    struct hash_uvec{
    size_t operator()(const arma::uvec & A)const{
        int b     = 378551;
        int a     = 63689;
        long hash = 0;
        for(int i = 0; i < A.n_rows; i++)
        {
            hash = hash * a + A(i);
            a    = a * b;
        }
        return hash;
    }
    };
    struct compare_uvec{
    bool operator()(const arma::uvec& __x,const arma::uvec& __y) const{
        if(__x.size()!=__y.size())return false;
        arma::uvec indices = arma::find( __x == __y );
        return indices.size() == __x.size() ;
    }
    };
    typedef __gnu_cxx::hash_map<arma::uvec,arma::uword,hash_uvec,compare_uvec> uvec_hash_map;
    RPEAC();
    void convert_to_atomic(
            const arma::umat&E,
            const arma::umat&P,
            arma::umat& aE,
            arma::umat& aP
            );

    void compute(
            const arma::umat& E,
            const arma::umat& P,
            arma::uvec& y
            );
    void show_atomic(arma::uvec& y);
    void convert_from_atomic(
            const arma::uvec& ay,
            arma::uvec& y
            );
private:
uvec_hash_map sig_hash_;
std::vector<arma::uword> index_map_;
};
}
#endif // RPEAC_H
