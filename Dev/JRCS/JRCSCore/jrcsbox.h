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
    virtual void initx(
            const MatPtr& xv,
            const MatPtr& xn,
            const CMatPtr& xc
            );
    virtual void reset_objw(const std::vector<float>&);
    static void set_boxes(std::vector<Cube::PtrLst>& cube_ptrlsts);
protected:
    void init_from_boxes();
    arma::fvec obj_prob_from_boxes(const Cube::PtrLst&, const MatPtr &vv);
    virtual void update_color_label();
private:
    static std::vector<Cube::PtrLst> cube_ptrlsts_;
    std::vector<std::shared_ptr<arma::gmm_diag>> obj_color_mode_;


};
}
#endif // JRCSBOX_H
