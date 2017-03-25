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
//    virtual void compute(void);
protected:
    virtual void prepare_compute();
    virtual void step_a(int i);
    virtual void step_b(void);
    virtual void calc_weighted(
            const arma::fmat&vv,
            const arma::fmat&vn,
            arma::Mat<uint8_t>&vc,
            const arma::mat& alpha,
            arma::fmat&wv,
            arma::fmat&wn,
            arma::Mat<uint8_t>&wc
            );
    virtual void updateRTforObj(
            const arma::uword start,
            const arma::uword end,
            arma::fmat& vv,
            arma::fmat& vn,
            arma::fmat& objv,
            arma::fmat& objn,
            arma::frowvec& alpha_colsum,
            arma::fmat& R,
            arma::fvec& t
            );

protected:
    virtual void update_color_label();
    virtual void reset_alpha();
    void init_from_boxes();
    arma::fvec obj_prob_from_boxes(const Cube::PtrLst&,const MatPtr &vv);
    void init_color_gmm(const Cube::PtrLst&,const MatPtr&,const CMatPtr&,GMMPtrLst&);
    void init_obj_prob(const Cube::PtrLst&,const MatPtr&,DMatPtr&);
    void init_color_prob(const CMatPtr&,DMatPtr&);
private:
    static std::vector<Cube::PtrLst> cube_ptrlsts_;
    GMMPtrLst color_gmm_lsts_;
    DMatPtrLst color_prob_lsts_;
    DMatPtrLst inbox_prob_lsts_;
    void debug_inbox_prob();
    void debug_color_prob();
};
}
#endif // JRCSBOX_H
