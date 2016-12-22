#ifndef SJRCSBASE_H
#define SJRCSBASE_H
#include "jrcscore_global.h"
#include <QCoreApplication>
#include "jrcsbase.h"
namespace JRCS{
//spectral joint registration and co-segmentation
class JRCSCORESHARED_EXPORT SJRCSBase:public JRCSBase
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
    SJRCSBase();
    virtual ~SJRCSBase(){}
    virtual std::string name()const{return "SJRCSBase";}
    virtual bool configure(Config::Ptr);
    virtual bool input_extra(
            const MeshBundle<DefaultMesh>::PtrList& inputs
            );
    virtual void initx(
            const MatPtr& xv,
            const MatPtr& xn,
            const CMatPtr& xc
            );//randomly initialize the X
    //Step A:
    //Calculate Alpha(Expectation)
    //Construct Residue Function
    //Step B:
    //Median of Residue Function
    //Step C:
    //Calculate Correlation Function
    //Alter Alpha
    //Update RT
    //Step D:
    //Update X
    //Update var pk
    virtual void compute(void)
    {
        std::cerr<<"preparing"<<std::endl;
        prepare_compute();
        while(!isEnd())
        {
            std::cerr<<"step a"<<std::endl;
            #pragma omp parallel for
            for( int i=0 ; i < vvs_ptrlst_.size() ; ++i )
            {
                step_a(i);
            }
            std::cerr<<"step b"<<std::endl;
            step_b();
            std::cerr<<"step c"<<std::endl;
            #pragma omp parallel for
            for( int i=0 ; i < vvs_ptrlst_.size() ; ++i )
            {
                step_c(i);
            }
            std::cerr<<"step d"<<std::endl;
            step_d();
            finish_steps();
        }
    }
protected:
    virtual void prepare_compute(void);
    virtual bool isEnd(void);
    virtual void step_a(int i);
    virtual void step_b(void);
    virtual void step_c(int i);
    virtual void step_d(void);
    virtual void finish_steps(void);
protected:
    virtual int evaluate_k();//evalute a proper number for x
    virtual void rand_sphere(
            arma::fmat& ov
            );//randomly sample point on sphere
    virtual void reset_prob();
    virtual void update_color_label();
    void prepare_for_residue_correlation(void);
    void truncate_alpha(const arma::mat&,arma::mat&,arma::vec&,arma::mat&,arma::vec&);
    void to_frame_func(const arma::mat&,const arma::vec&,const arma::vec&,arma::vec&);
    void to_vox_func(int index,const arma::vec& frame_func,arma::vec& vox_func);
    void proj_and_rebuild(int index,const arma::vec& vox_func,arma::vec& re_vox_func);
    void to_pix_func(int index,const arma::vec& vox_func,arma::vec& pix_func);
    void to_model_func(const arma::mat&,const arma::vec&,const arma::vec&,arma::vec&);
    double res_energy(const arma::vec& func0,const arma::vec& func1);
    void max_cir_cor(const arma::vec&,const arma::vec&,arma::uword& max_cir);//maximum of circular correlation
    void cir_mat(const arma::uword& cir,arma::mat& alpha);//perform cirulation on alpha
    void calc_weighted(
            const arma::fmat&vv,
            const arma::fmat&vn,
            arma::Mat<uint8_t>&vc,
            const arma::mat& alpha,
            arma::fmat&wv,
            arma::fmat&wn,
            arma::Mat<uint8_t>&wc
            );
    void reset_x(void);
    void reset_var_p(void);
    virtual void calc_obj(void);
protected:
    //iteration count
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

    //result correspondece
    using JRCSBase::alpha_ptrlst_;

    //use residue
    bool use_res_;
    arma::uword res_step_;
    arma::uword rebuild_k_;
    arma::uword base_k_;
    arma::uword res_act_freq_;
    //circulated indicating functions
    arma::mat funcs_;
    //circulating position
    arma::uvec cir_pos_;
    //circulated indices
    arma::umat cir_indices_;
    //cir frame
    arma::uword cir_frame_;
    //cir value
    arma::uvec cir_value_;

    //residue function
    arma::mat res_;
    arma::vec median_res_;
    MeshList inputs_;
    DMatPtrLst bases_;

    //weighted V
    using JRCSBase::wvs_ptrlst_;
    using JRCSBase::wns_ptrlst_;
    using JRCSBase::wcs_ptrlst_;

    //latent model centroid
    using JRCSBase::xv_ptr_;
    using JRCSBase::xn_ptr_;
    using JRCSBase::xc_ptr_;

    //transformed centroid
    MatPtrLst  xtv_ptrlst_;

    arma::mat vvar_;
    //latent model parameter
    using JRCSBase::x_p_;      //probability
    using JRCSBase::x_invvar_; //inversed variance
    using JRCSBase::beta_;      //noise portion

    //pre-defined object info
    using JRCSBase::obj_num_;
    std::vector<arma::uword> obj_range_;//start and end index for objects in centroid model
    using JRCSBase::obj_pos_;

    //minimization object
    double obj_;
    //configuration
    Config::Ptr config_;
    //verbose level
    using JRCSBase::verbose_;
};
}
#endif // SJRCSBASE_H
