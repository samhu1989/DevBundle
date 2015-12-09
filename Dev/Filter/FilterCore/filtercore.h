#ifndef FILTERCORE_H
#define FILTERCORE_H
#include "common.h"
#include "filtercore_global.h"
namespace Filter {
template<typename M>
class FilterBase
{
public:
    inline void extract(M&m)
    {
        if(!initCompute())return;
        applyFilter(m);
        deinitCompute();
    }
protected:
    virtual bool initCompute();
    virtual void deinitCompute();
    virtual void applyFilter(M&m)=0;
};
}
#include "filtercore.hpp"
#endif // FILTERCORE_H
