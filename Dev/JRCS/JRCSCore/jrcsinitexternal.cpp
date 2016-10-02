#include "jrcsinitexternal.h"
namespace JRCS{
bool JRCSInitExternal::configure(Config::Ptr config)
{
    return true;
}
bool JRCSInitExternal::init_with_label(
        const int k,
        const MatPtrLst& vv,
        const MatPtrLst& vn,
        const CMatPtrLst& vc,
        const LCMatPtrLst& vlc,
        const LMatPtrLst& vl,
        int verbose
        )
{
    if(verbose)std::cerr<<"Init by "<<name()<<std::endl;
    return true;
}
}
