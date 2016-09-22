#ifndef JRCSINITOUTSIDE_H
#define JRCSINITOUTSIDE_H
#include "jrcsinitbase.h"
namespace JRCS {
class JRCSCORESHARED_EXPORT JRCSInitExternal:public JRCSInitBase
{
public:
    JRCSInitExternal(MatPtrLst& a,arma::fvec& p):external_alpha_(a),external_obj_prob_(p),JRCSInitBase(){}
    virtual bool configure(Config::Ptr config);
    virtual bool init_with_label(
            const int k,
            const MatPtrLst& vv,
            const MatPtrLst& vn,
            const CMatPtrLst& vc,
            const LCMatPtrLst& vlc,
            const LMatPtrLst& vl,
            int verbose
            );
    virtual void getAlpha(MatPtrLst& alpha){alpha = external_alpha_;}
    virtual void getObjProb(arma::fvec& obj_prob){obj_prob = external_obj_prob_;}
protected:
private:
    MatPtrLst& external_alpha_;
    arma::fvec& external_obj_prob_;
};
}
#endif // JRCSINITOUTSIDE_H
