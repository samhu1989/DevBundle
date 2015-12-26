#ifndef BLOCKBASEDFEATURE_H
#define BLOCKBASEDFEATURE_H
#include <armadillo>
namespace Feature{
template<typename Mesh>
class BlockBasedFeature
{
public:
    BlockBasedFeature();
    void extract(const Mesh&,arma::vec&);
};
}
#endif // BLOCKBASEDFEATURE_H
