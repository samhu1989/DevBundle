#ifndef JRCSCORE_H
#define JRCSCORE_H
#include "jrcscore_global.h"
#include <armadillo>
#include <memory>
#include <QCoreApplication>
#include "jrcsinitbase.h"
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
    typedef std::shared_ptr<arma::uvec> LMatPtr;
    typedef std::vector<LMatPtr> LMatPtrLst;
    typedef struct{
        float R[9];
        float t[3];
    }T;
    typedef std::vector<T> Ts;
    typedef std::vector<Ts> TsLst;
    typedef enum{
        ObjOnly,
        ObjPointDist
    }CompatibilityType;
    typedef enum{
        All,
        Alpha,
        Beta,
        Gamma
    }RotationType;
    JRCSBase(){arma::arma_rng::set_seed_random();}
    virtual ~JRCSBase(){}

    virtual bool configure(Config::Ptr config);

    virtual inline void enable_smooth(bool enable=true){smooth_enabled_=enable;}
    virtual inline void set_smooth_weight(float w){smooth_w_=w;}
    virtual inline void set_max_smooth_iter(int max_iter){max_smooth_iter_=max_iter;}
    virtual inline void set_max_iter(int max_iter){max_iter_ = max_iter;}
    virtual inline void set_max_init_iter(int max){max_init_iter_=max;}
    virtual inline void set_debug_path(const std::string& path){debug_path_=path;}
    virtual inline void set_mu_type(const CompatibilityType& type){mu_type_=type;}
    virtual inline void set_rt_type(const RotationType& type){rttype_=type;}
    virtual inline int  get_iter_num(void){return iter_count_;}
    virtual inline int  get_max_iter(void){return max_iter_;}
    virtual inline void get_rt(TsLst& rt){rt = rt_lst_;}
    virtual inline int  get_obj_num(void){return obj_num_;}

    virtual void get_label(std::vector<arma::uvec>&);
    virtual void input(
            const MatPtrLst& vv,
            const MatPtrLst& vn,
            const CMatPtrLst& vc,
            const LCMatPtrLst& vl
            );
    virtual void input_with_label(
            const MatPtrLst& vv,
            const MatPtrLst& vn,
            const CMatPtrLst& vc,
            const LCMatPtrLst& vlc,
            const LMatPtrLst& vl
            );
    virtual void resetw(
            const MatPtrLst& wv,
            const MatPtrLst& wn,
            const CMatPtrLst& wc
            );
    virtual inline void reset_n( int n = 4 ){
        obj_num_ = n + 1;
        if( obj_num_ < 2 )throw std::logic_error("works for at least two object number");
    }
    virtual void reset_objw(const std::vector<float>&);
    virtual int evaluate_k();//evalute a proper x
    virtual void initx(
            const MatPtr& xv,
            const MatPtr& xn,
            const CMatPtr& xc
            );//randomly initialize the X

    virtual void reset_rt();
    virtual void reset_iteration();
    virtual void compute()
    {
        reset_alpha();
        reset_prob();
        while(!isEnd())
        {
            computeOnce();
            update_color_label();
            QCoreApplication::processEvents();
            ++iter_count_;
        }
    }
protected:
    virtual void obj_only(arma::mat&mu);
    virtual void obj_point_dist(arma::mat&mu);
    virtual void computeCompatibility(arma::mat& mu);
    virtual void computeOnce();
    virtual bool isEnd();
    virtual void reset_obj_vn(
            float radius,
            arma::fvec& pos,
            arma::fmat& ov,
            arma::fmat& on
            );
    virtual void reset_obj_c(
            arma::Mat<uint8_t>& oc
            );
    virtual void reset_alpha();
    virtual void reset_prob();
    static void rand_sphere(
            arma::fmat& ov
            );
    virtual void update_color_label();
protected:
    //object number
    int obj_num_;

    //termination criteria
    int iter_count_;

    //max radius
    float max_obj_radius_;

    //max iteration number
    int max_iter_;

    //max init iteration number
    int max_init_iter_;

    //smooth iter number
    int max_smooth_iter_;

    //input observation
    MatPtrLst vvs_ptrlst_;
    MatPtrLst vns_ptrlst_;
    CMatPtrLst vcs_ptrlst_;
    LCMatPtrLst vls_ptrlst_;
    LMatPtrLst vll_ptrlst_;
    TsLst rt_lst_;

    //results
    MatPtrLst alpha_ptrlst_;

    //weighted V
    MatPtrLst  wvs_ptrlst_;
    MatPtrLst  wns_ptrlst_;
    CMatPtrLst wcs_ptrlst_;

    //latent model centroid
    MatPtr xv_ptr_;
    MatPtr xn_ptr_;
    CMatPtr xc_ptr_;

    //transformed latent model centroid
    arma::fmat xtv_;
    arma::fmat xtn_;
    arma::Mat<uint8_t> xtc_;

    //pre-defined class label for each centroid
    arma::uvec  obj_label_;
    arma::fvec  obj_prob_;
    arma::fmat  obj_pos_;

    MatPtrLst   objv_ptrlst_;
    MatPtrLst   objn_ptrlst_;
    CMatPtrLst  objc_ptrlst_;

    //latent model parameter
    arma::frowvec x_p_;
    arma::frowvec x_invvar_;

    //x compatiblity
    CompatibilityType mu_type_;
    arma::mat mu_;

    //sum of latent model
    arma::fmat xv_sum_;
    arma::fmat xn_sum_;
    arma::fmat xc_sum_;

    arma::frowvec var_sum;
    arma::frowvec alpha_sum;
    arma::frowvec alpha_sumij;

    int verbose_;
    std::string debug_path_;
    bool smooth_enabled_;
    RotationType rttype_;
    float beta_;
    float smooth_w_;
    bool init_alpha_;
private:
    std::shared_ptr<JRCSInitBase> init_;
    Config::Ptr config_;
};
}
#endif // JRCSCORE_H
