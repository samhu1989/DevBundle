#include "jrcscore.h"
namespace JRCS{

void JRCSBase::reset_iteration()
{
    iter_count_ = 0;
    obj_num_ = 3;
}

void JRCSBase::input(
        const MatPtrLst& vv,
        const MatPtrLst& vn,
        const CMatPtrLst& vc,
        const LCMatPtrLst& vl
        )
{
    vvs_ptrlst_ = vv;
    vns_ptrlst_ = vn;
    vcs_ptrlst_ = vc;
    vls_ptrlst_ = vl;
}

void JRCSBase::resetw(
        const MatPtrLst& wv,
        const MatPtrLst& wn,
        const CMatPtrLst& wc
        )
{
    wvs_ptrlst_ = wv;
    wns_ptrlst_ = wn;
    wcs_ptrlst_ = wc;
}

void JRCSBase::init_x()
{

}

int JRCSBase::evaluate_k()
{

}

void JRCSBase::computeOnce()
{
    stepE();
    stepM();
}

void JRCSBase::stepE()
{
    ;
}

void JRCSBase::stepM()
{
    ;
}

bool JRCSBase::isEnd()
{
    ;
}
}
