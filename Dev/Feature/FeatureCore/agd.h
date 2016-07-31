#ifndef AGD_H
#define AGD_H
#include <armadillo>
#include "common.h"
namespace Feature{
template<typename Mesh>
class AGD
{
public:
    class Pair{
    public:
        Pair(uint16_t f,uint16_t s,double d):first_(f),second_(s),dist_(d){}
        uint16_t first_;
        uint16_t second_;
        double dist_;
        bool operator < (const Pair &a) const {
            return dist_ > a.dist_;//little first
        }
    };
    void extract(const Mesh&,arma::vec&);
    void extract(const VoxelGraph<Mesh>&,arma::vec&);
    void extract_hash(const VoxelGraph<Mesh>&,arma::vec&);
    void extract_slow(const VoxelGraph<Mesh>&,arma::vec&);
};
}
#endif // AGD_H
