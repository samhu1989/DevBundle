#ifndef JRCSBOX_H
#define JRCSBOX_H
#include "jrcscore_global.h"
#include "jrcsbilateral.h"
#include "cube.h"
namespace JRCS{
class JRCSCORESHARED_EXPORT JRCSBox:public JRCSBilateral
{
public:
    using JRCSBase::DMatPtr;
    using JRCSBase::CMatPtrLst;
    typedef Common::Cube Cube;
    typedef std::shared_ptr<arma::gmm_diag> GMMPtr;
    typedef std::vector<GMMPtr> GMMPtrLst;

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
    virtual void update_color_label();
    virtual void reset_alpha();
    void init_from_boxes();
    arma::fvec obj_prob_from_boxes(const Cube::PtrLst&, const MatPtr &vv);
    void init_color_gmm(const Cube::PtrLst&, const MatPtr&, const CMatPtr&, GMMPtrLst&);
    void init_obj_prob(const Cube::PtrLst&,const MatPtr&,DMatPtrLst&);
private:
    static std::vector<Cube::PtrLst> cube_ptrlsts_;
    std::vector<GMMPtrLst> color_gmm_lsts_;
    std::vector<DMatPtrLst> obj_prob_lsts_;
};
}
#endif // JRCSBOX_H
