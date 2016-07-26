#ifndef AGD_H
#define AGD_H
#include <armadillo>
#include "common.h"
namespace Feature{
template<typename Mesh>
class AGD
{
public:
    void extract(const Mesh&,arma::vec&);
    void extract(const VoxelGraph<Mesh>&,arma::vec&);
};
}
#endif // AGD_H
