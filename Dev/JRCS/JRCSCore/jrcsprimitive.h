#ifndef JRCSPRIMITIVE_H
#define JRCSPRIMITIVE_H
#include "jrcsbilateral.h"
#include <ext/hash_map>"
namespace JRCS{
struct  Plate{
    typedef std::shared_ptr<Plate> Ptr;
private:
    arma::fmat R_;
    arma::fvec t_;
    //voting space
    arma::fcube param_z_u_;
    arma::fcube param_x_f_;
    arma::fcube param_x_b_;
    arma::fcube param_y_f_;
    arma::fcube param_y_b_;
    //color space
    arma::fvec  rgb_z_u_;
    arma::fvec  rgb_x_f_;
    arma::fvec  rgb_x_b_;
    arma::fvec  rgb_y_f_;
    arma::fvec  rgb_y_b_;
    //parameter
    arma::fvec  res_z_u_;
    arma::fvec  res_x_f_;
    arma::fvec  res_x_b_;
    arma::fvec  res_y_f_;
    arma::fvec  res_y_b_;
};
class JRCSCORESHARED_EXPORT JRCSPrimitive:public JRCSBilateral
{
public:
    JRCSPrimitive();
    virtual ~JRCSPrimitive(){}
    virtual std::string name()const{return "JRCSPrimitive";}
    virtual void initx(
            const MatPtr& xv,
            const MatPtr& xn,
            const CMatPtr& xc
            );//randomly initialize the X
protected:
    virtual void reset_obj_vn(
            float radius,
            arma::fvec& pos,
            arma::fmat& ov,
            arma::fmat& on
            );
    virtual void compute(void);
protected:
    virtual void prepare_primitive();
    virtual void finish_primitive();
    //calculate alpha
    //update r t
    //voting
    void step_1(int i);
    //extracting new planes
    //updating var and p
    void step_2(void);
    virtual bool isEnd_primitive(void);
private:
    std::vector<Plate::Ptr> plate_ptrlst_;
};
}
#endif // JRCSHOUGH_H
