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
    virtual void prepare_alpha_operation(int i);
    virtual void alpha_operation(int i);
    virtual void alpha_operation_a(int i);
    virtual void alpha_operation_b(int i);
    virtual void rand_sphere(
            arma::fmat& ov
            );
private:
    std::vector<arma::mat> init_alpha_value;
};
}

#endif // JRCSAOPT_H
