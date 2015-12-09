#ifndef FILTERCORE_HPP
#define FILTERCORE_HPP
#include "filtercore.h"
namespace Filter {
template<typename M>
bool FilterBase<M>::initCompute()
{
    return false;
}
template<typename M>
void FilterBase<M>::deinitCompute()
{

}
}
#endif // FILTERCORE_HPP

