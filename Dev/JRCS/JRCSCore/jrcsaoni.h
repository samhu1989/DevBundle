#ifndef JRCSAONI_H
#define JRCSAONI_H
#include "jrcsbase.h"
//JRCS average on neareast input
namespace JRCS {
class JRCSCORESHARED_EXPORT JRCSAONI : public JRCSBase
{
public:
    JRCSAONI();
    virtual std::string name()const{return "JRCSAONI";}
protected:
    virtual void computeOnce();
    virtual void alpha_operation(int i){}
    virtual void prepare_alpha_operation(int i){}
};
}
#endif // JRCSAONI_H
