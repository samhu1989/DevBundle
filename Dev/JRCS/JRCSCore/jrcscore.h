#ifndef JRCSCORE_H
#define JRCSCORE_H
#include "jrcscore_global.h"
#include <armadillo>
#include <memory>
namespace JRCS
{
class JRCSCORESHARED_EXPORT JRCSBase
{
public:
    typedef std::shared_ptr<arma::fmat> MatPtr;
    typedef std::vector<MatPtr> MatPtrLst;
    typedef std::shared_ptr<arma::Mat<uint8_t>> CMatPtr;
    typedef std::vector<CMatPtr> CMatPtrLst;
    typedef std::shared_ptr<arma::Col<uint32_t>> LCMatPtr;
    typedef std::vector<LCMatPtr> LCMatPtrLst;
    JRCSBase(){}
    virtual ~JRCSBase(){}
    virtual void input(
            const MatPtrLst& vv,
            const MatPtrLst& vn,
            const CMatPtrLst& vc,
            const LCMatPtrLst& vl
            );
    virtual void resetw(
            const MatPtrLst& wv,
            const MatPtrLst& wn,
            const CMatPtrLst& wc
            );
    virtual inline void reset_n(int n=4){obj_num_=n;}
    virtual int evaluate_k();//evalute a proper x
    virtual void init_x();//randomly initialize the X
    virtual void reset_iteration();
    virtual void compute()
    {
        while(!isEnd())
        {
            computeOnce();
            ++iter_count_;
        }
    }
protected:
    virtual void computeOnce();
    virtual void stepE();
    virtual void stepM();
    virtual bool isEnd();
protected:
    //object number
    int obj_num_;

    //termination criteria
    int iter_count_;

    //input observation
    MatPtrLst vvs_ptrlst_;
    MatPtrLst vns_ptrlst_;
    CMatPtrLst vcs_ptrlst_;
    LCMatPtrLst vls_ptrlst_;

    //results
    MatPtrLst alpha_ptrlst_;

    //weighted V
    MatPtrLst  wvs_ptrlst_;
    MatPtrLst  wns_ptrlst_;
    CMatPtrLst wcs_ptrlst_;

    //latent model centroid
    MatPtr x_ptr_;
    MatPtr xn_ptr_;
    CMatPtr xc_ptr_;
    arma::uvec obj_label_;
    //pre-defined class label for each centroid

    //latent model parameter
    arma::frowvec x_p_;
    arma::frowvec x_invvar_;

    arma::fvec obj_p_;

    MatPtrLst obj_ptrlst_;
    MatPtrLst objn_ptrlst_;
    //probability
};
}
#endif // JRCSCORE_H
