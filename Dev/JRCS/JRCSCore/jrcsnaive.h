#ifndef JRCSNAIVE_H
#define JRCSNAIVE_H
#include "jrcscore_global.h"
#include "jrcsbase.h"
namespace JRCS
{
class JRCSCORESHARED_EXPORT JRCSNaive:public JRCSBase
{
public:
    JRCSNaive():JRCSBase(){}
protected:
    void computeOnce();
};
}
#endif // JRCSNAIVE_H
