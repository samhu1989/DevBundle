#ifndef JRCSBOX_H
#define JRCSBOX_H
#include "jrcscore_global.h"
#include "jrcsbilateral.h"
#include "cube.h"
namespace JRCS{
class JRCSCORESHARED_EXPORT JRCSBox:public JRCSBilateral
{
public:
    typedef Common::Cube Cube;
    JRCSBox();
    virtual ~JRCSBox(){}
    virtual std::string name()const{return "JRCSBox";}
    static void set_boxes(std::vector<Cube::PtrLst>& cube_ptrlsts);
    virtual void initx(
            const MatPtr& xv,
            const MatPtr& xn,
            const CMatPtr& xc
            );//initialize the X
protected:
    virtual void compute(void);
    virtual void reset_obj_vn(
            float radius,
            arma::fvec& pos,
            arma::fmat& ov,
            arma::fmat& on
            );
protected:
    void update_objective();
private:
    static std::vector<Cube::PtrLst> cube_ptrlsts_;
};
}
#endif // JRCSBOX_H
