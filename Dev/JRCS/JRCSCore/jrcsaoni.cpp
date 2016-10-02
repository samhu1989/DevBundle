#include "jrcsaoni.h"
namespace JRCS {
JRCSAONI::JRCSAONI():JRCSBase()
{

}
void JRCSAONI::computeOnce()
{
    //reset sum
    xv_sum_.fill(0.0);
    xn_sum_.fill(0.0);
    xc_sum_.fill(0.0);
    var_sum.fill(0.0);
    alpha_sum.fill(0.0);
    alpha_sumij.fill(0.0);

    //reset transformed latent center
    xtc_ = *xc_ptr_;

    for( int i = 0 ; i < vvs_ptrlst_.size() ; ++ i )
    {
    //input vertices normal color
        arma::fmat& vv_ = *vvs_ptrlst_[i];
        arma::fmat& vn_ = *vns_ptrlst_[i];
        arma::Mat<uint8_t>& vc_ = *vcs_ptrlst_[i];
    //update alpha
        arma::fmat& alpha = *alpha_ptrlst_[i];
        Ts& rt = rt_lst_[i];
    //transform the object
        #pragma omp parallel for
        for(int o = 0 ; o < obj_num_ ; ++o )
        {
            arma::fmat R(rt[o].R,3,3,false,true);
            arma::fvec t(rt[o].t,3,false,true);
            arma::fmat& objv = *objv_ptrlst_[o];
            arma::fmat& objn = *objn_ptrlst_[o];
            objv = R*objv;
            objv.each_col() += t;
            objn = R*objn;
        }
    //update RT
    }
    //update X

    //fix the X center position

}
}
