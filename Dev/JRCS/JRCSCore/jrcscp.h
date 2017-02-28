#ifndef JRCSCP_H
#define JRCSCP_H
/*
 * JRCS Primitive with more constraints
 */
#include "jrcsprimitive.h"
namespace JRCS {
struct JRCSCORESHARED_EXPORT CPlate:public Plate
{
    typedef enum{
        H, // horizontal
        V  // vertical
    }T;
    CPlate();
    CPlate(const arma::fmat& v,
        const arma::fmat& n,
        const arma::Mat<uint8_t>& c,
        const arma::fvec &pos
    );
    virtual void local_translate(
        const arma::fvec& t,
        Plate& result
    );
    virtual void scale(
        const arma::fvec &s,
        Plate &result
    );
private:
    T type_;
};
class JRCSCORESHARED_EXPORT JRCSCP:public JRCSPrimitive
{
public:
    JRCSCP();
    virtual ~JRCSCP(){}
    virtual std::string name()const{return "JRCSCP";}
    virtual void initx(
            const MatPtr& xv,
            const MatPtr& xn,
            const CMatPtr& xc
            );//randomly initialize the X
protected:
    virtual void prepare_primitive();
private:
    using JRCSPrimitive::plate_ptrlst_;
    using JRCSPrimitive::plate_t_ptrlst_;
    using JRCSPrimitive::plate_num_for_obj_;
    using JRCSPrimitive::point_num_for_plate_;
};
}
#endif // JRCSCP_H
