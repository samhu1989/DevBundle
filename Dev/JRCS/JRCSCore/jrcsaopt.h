#ifndef JRCSAOPT_H
#define JRCSAOPT_H
#include "jrcsaoni.h"
//JRCS with operation on alpha
namespace JRCS {
class JRCSCORESHARED_EXPORT JRCSAOPT: public JRCSAONI
{
public:
    JRCSAOPT();
    virtual std::string name()const{return "JRCSAOPT";}
protected:
    virtual void alpha_operation(int i);
    virtual void rand_sphere(
            arma::fmat& ov
            );
};
}

#endif // JRCSAOPT_H
