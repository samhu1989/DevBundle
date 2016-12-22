#ifndef JRCSBILATERAL_H
#define JRCSBILATERAL_H
#include "sjrcsbase.h"
#include "jrcsbilateral.h"
namespace JRCS{
class JRCSCORESHARED_EXPORT JRCSBilateral:public SJRCSBase
{
public:
    using JRCSBase::MatPtr;
    using JRCSBase::MatPtrLst;
    using JRCSBase::DMatPtr;
    using JRCSBase::DMatPtrLst;
    using JRCSBase::CMatPtr;
    using JRCSBase::CMatPtrLst;
    using JRCSBase::LCMatPtr;
    using JRCSBase::LCMatPtrLst;
    using JRCSBase::LMatPtr;
    using JRCSBase::LMatPtrLst;
    typedef std::shared_ptr<arma::vec> VecPtr;
    typedef std::vector<VecPtr> VecPtrLst;
    using JRCSBase::T;
    using JRCSBase::Ts;
    using JRCSBase::TsLst;
    using JRCSBase::RotationType;
    typedef std::vector<MeshBundle<DefaultMesh>::Ptr> MeshList;
public:
    JRCSBilateral();
    virtual ~JRCSBilateral(){}
    virtual std::string name()const{return "JRCSBilateral";}
    virtual bool configure(Config::Ptr);
    virtual void compute(void);
    virtual bool input_extra(const MeshBundle<DefaultMesh>::PtrList& inputs);
protected:
    virtual void prepare_compute(void);
    virtual void step_a(int i);
    virtual void step_b(void);
    using SJRCSBase::finish_steps;
protected:
    virtual void reset_prob();
    virtual void calc_obj();
    void calc_weighted(
            const arma::fmat&vv,
            const arma::fmat&vn,
            const arma::fmat&vf,
            arma::Mat<uint8_t>&vc,
            const arma::mat& alpha,
            arma::fmat&wv,
            arma::fmat&wn,
            arma::fmat&wf,
            arma::Mat<uint8_t>&wc
            );
    void rescale_feature();
//    void proj_and_rebuild(const int i,const arma::vec& vox_func,arma::vec& re_vox_func);
//    void smooth_on_alpha(const int i, arma::mat& alpha);
protected:
    using JRCSBase::iter_count_;

    //max iteration number
    using JRCSBase::max_iter_;

    //type of rotation
    using JRCSBase::rttype_;

    //input observation
    using JRCSBase::vvs_ptrlst_;
    using JRCSBase::vns_ptrlst_;
    using JRCSBase::vcs_ptrlst_;
    using JRCSBase::vls_ptrlst_; //label color
    using JRCSBase::vll_ptrlst_; //labels
    MatPtrLst vfs_ptrlst_;

    //
//    using SJRCSBase::bases_;
//    using SJRCSBase::base_k_;
//    using SJRCSBase::rebuild_k_;
    //result correspondece
    using JRCSBase::alpha_ptrlst_;

    //weighted V
    using JRCSBase::wvs_ptrlst_;
    using JRCSBase::wns_ptrlst_;
    using JRCSBase::wcs_ptrlst_;
    MatPtrLst wfs_ptrlst_;

    //latent model centroid
    using JRCSBase::xv_ptr_;
    using JRCSBase::xn_ptr_;
    using JRCSBase::xc_ptr_;
    MatPtr xf_ptr_;

    //transformed centroid
    using SJRCSBase::xtv_ptrlst_;
    MatPtrLst xtn_ptrlst_;

    using SJRCSBase::vvar_;
    arma::mat fvar_;
    //latent model parameter
    using SJRCSBase::x_p_;      //vertex probability
    arma::rowvec xv_invvar_; //vertex inversed variance
    arma::rowvec xf_invvar_; //normal inversed variance
    using JRCSBase::beta_;      //noise portion

    //pre-defined object info
    arma::uword f_dim_ ; //feature dimenssion
    using JRCSBase::obj_num_;
    using SJRCSBase::obj_range_;//start and end index for objects in centroid model
    using JRCSBase::obj_pos_;

    //minimization object
    using SJRCSBase::obj_;
    double obj_f;
    double obj_v;
    //configuration
    Config::Ptr config_;
    //verbose level
    using JRCSBase::verbose_;
};
}
#endif // JRCSBILATERAL_H
