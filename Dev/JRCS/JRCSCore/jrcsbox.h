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
    virtual void reset_rt();
    virtual void reset_objw(const std::vector<float>&);
    static void set_boxes(std::vector<Cube::PtrLst>& cube_ptrlsts);
    static void add_boxes(std::vector<Cube::PtrLst>& cube_ptrlsts);
    virtual void get_label(std::vector<arma::uvec>&);
    virtual void get_order(std::vector<arma::uvec>&);
//    virtual void compute(void);
protected:
    virtual void prepare_compute();
    virtual void step_a(int i);
    virtual void finish_steps();

protected:
    void update_from_cube(void);
    virtual void update_color_label();
    virtual void reset_alpha();
    void init_from_boxes();
    arma::fvec obj_prob_from_boxes(const Cube::PtrLst&,const MatPtr &vv);
    void init_color_gmm(const Cube::PtrLst&,const MatPtr&,const CMatPtr&,GMMPtrLst&);
    void init_obj_prob(const Cube::PtrLst&,const MatPtr&,DMatPtr&);
    void init_color_prob(const CMatPtr&,DMatPtr&);
    void init_obj_x(
            int obj_idx,
            arma::fmat& objv,
            arma::fmat& objn,
            arma::Mat<uint8_t>& objc,
            const arma::fvec& pos
            );
private:
    static std::vector<Cube::PtrLst> cube_ptrlsts_;
    static bool update_cube_;
    arma::uword cube_init_frame_;
    std::vector<std::vector<arma::uword>> obj_cube_index;
    GMMPtrLst color_gmm_lsts_;
    DMatPtrLst color_prob_lsts_;
    DMatPtrLst inbox_prob_lsts_;
    void debug_inbox_prob();
    void debug_color_prob();
    void debug_alpha(int i);
};
}
#endif // JRCSBOX_H
