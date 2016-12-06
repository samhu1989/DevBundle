#include "sjrcsbase.h"
namespace JRCS{
SJRCSBase::SJRCSBase()
{
    arma::arma_rng::set_seed_random();
}
bool SJRCSBase::configure(Config::Ptr config)
{
    config_ = config;
}
void SJRCSBase::input(
        const MatPtrLst& vv,
        const MatPtrLst& vn,
        const CMatPtrLst& vc,
        const LCMatPtrLst& vl
        )
{
    ;
}
void SJRCSBase::resetw(
        const MatPtrLst& wv,
        const MatPtrLst& wn,
        const CMatPtrLst& wc
        )
{
    ;
}
void SJRCSBase::initx(
        const MatPtr& xv,
        const MatPtr& xn,
        const CMatPtr& xc
        )
{
    ;
}
bool SJRCSBase::isEnd(void)
{
    ;
}
void SJRCSBase::prepare_compute(void)
{
    ;
}
void SJRCSBase::parallel_compute(int i)
{
    ;
}
void SJRCSBase::finish_compute(void)
{
    ;
}
}
