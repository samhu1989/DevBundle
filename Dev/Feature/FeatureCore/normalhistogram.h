#ifndef NORMALHISTOGRAM_H
#define NORMALHISTOGRAM_H
#include "common.h"
namespace Feature
{
template<typename Mesh>
class NormalHistogram
{
public:
    NormalHistogram();
    void extract(Mesh&,arma::vec&);
protected:
};
}
#include "normalhistogram.hpp"
#endif // NORMALHISTOGRAM_H
