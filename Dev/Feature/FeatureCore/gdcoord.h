#ifndef GDCOORD_H
#define GDCOORD_H
#include "common.h"
namespace Feature{
//Geodesic Distance Coordinates
template<typename Mesh>
class GDCoord
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
    void extract(const VoxelGraph<Mesh>& graph, const arma::uvec& axis_coord ,arma::fmat&);
};
}
#endif // GDCOORD_H
