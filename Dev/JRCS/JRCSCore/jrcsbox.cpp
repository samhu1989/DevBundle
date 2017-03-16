#include "jrcsbox.h"
namespace JRCS{

std::vector<JRCSBox::Cube::PtrLst> JRCSBox::cube_ptrlsts_;

void JRCSBox::set_boxes(std::vector<Cube::PtrLst>& cube_ptrlsts)
{
    cube_ptrlsts_ = cube_ptrlsts;
}

JRCSBox::JRCSBox()
{
    ;
}

void JRCSBox::initx(
        const MatPtr& xv,
        const MatPtr& xn,
        const CMatPtr& xc
        )
{
    ;
}

void JRCSBox::compute(void)
{
    ;
}

void JRCSBox::reset_obj_vn(
        float radius,
        arma::fvec& pos,
        arma::fmat& ov,
        arma::fmat& on
        )
{
    ;
}

void JRCSBox::update_objective()
{
    ;
}

}
