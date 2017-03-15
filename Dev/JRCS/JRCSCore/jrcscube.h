#ifndef JRCSCUBE_H
#define JRCSCUBE_H
#include "jrcsbilateral.h"
#include "cube.h"
namespace JRCS {
class JRCSCORESHARED_EXPORT JRCSCube:public JRCSBilateral
{
public:
    typedef Common::Cube Cube;
    JRCSCube();
    virtual ~JRCSCube(){}
    virtual std::string name()const{return "JRCSCube";}
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
    virtual void reset_alpha_cube();
    virtual void reset_prob_cube();
    virtual void prepare_cube();
    virtual void finish_cube();
    void update_objective();
    void output_debug();
    //calculate alpha
    //update r t
    //voting
    virtual void step_1(int i);
    virtual void updateRTforObj(
            const int start,
            const int end,
            arma::frowvec& colsum,
            arma::fmat& R,
            arma::fvec& t,
            Cube::PtrLst cube_ptrlst
            );
    //extracting new planes
    //updating var and p
    virtual void step_2(void);
    virtual bool isEnd_cube(void);
protected:
    Cube::PtrLst cube_ptrlst_;
    std::vector<Cube::PtrLst> cube_t_ptrlst_;
};
}
#endif // JRCSCUBE_H
