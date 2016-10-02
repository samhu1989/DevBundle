#ifndef JRCSAONI_H
#define JRCSAONI_H
#include "jrcsbase.h"
namespace JRCS {
class JRCSCORESHARED_EXPORT JRCSAONI : public JRCSBase
{
public:
    JRCSAONI();
    virtual std::string name()const{return "JRCSAONI";}
protected:
    virtual void computeOnce();
};
}
#endif // JRCSAONI_H
