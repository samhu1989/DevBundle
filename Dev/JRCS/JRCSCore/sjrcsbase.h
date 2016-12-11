#ifndef SJRCSBASE_H
#define SJRCSBASE_H
#include "jrcscore_global.h"
#include "common.h"
#include <QCoreApplication>
namespace JRCS{
//spectral joint registration and co-segmentation
class JRCSCORESHARED_EXPORT SJRCSBase
{
public:
    typedef std::shared_ptr<arma::fmat> MatPtr;
    typedef std::vector<MatPtr> MatPtrLst;
    typedef std::shared_ptr<arma::mat> DMatPtr;
    typedef std::vector<DMatPtr> DMatPtrLst;
    typedef std::shared_ptr<arma::Mat<uint8_t>> CMatPtr;
    typedef std::vector<CMatPtr> CMatPtrLst;
    typedef std::shared_ptr<arma::Col<uint32_t>> LCMatPtr;
    typedef std::vector<LCMatPtr> LCMatPtrLst;
    typedef std::shared_ptr<arma::uvec> LMatPtr;
    typedef std::vector<LMatPtr> LMatPtrLst;
    typedef std::shared_ptr<arma::vec> VecPtr;
    typedef std::vector<VecPtr> VecPtrLst;
    typedef struct{
        float R[9];
        float t[3];
    }T;
    typedef std::vector<T> Ts;
    typedef std::vector<Ts> TsLst;
    typedef enum{
        All,
        Alpha,
        Beta,
        Gamma
    }RotationType;
public:
    SJRCSBase();
    virtual ~SJRCSBase(){}
    virtual std::string name()const{return "SJRCSBase";}
    virtual bool configure(Config::Ptr);
    //input observation from outside
    virtual void input(
            const MatPtrLst& vv,
            const MatPtrLst& vn,
            const CMatPtrLst& vc,
            const LCMatPtrLst& vl
            );
    //pass allocated from outside
    virtual void resetw(
            const MatPtrLst& wv,
            const MatPtrLst& wn,
            const CMatPtrLst& wc
            );
    //randomly initialize the X
    virtual void initx(
            const MatPtr& xv,
            const MatPtr& xn,
            const CMatPtr& xc
            );
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
        prepare_compute();
        while(!isEnd())
        {
            #pragma omp parallel for
            for( int i=0 ; i < vvs_ptrlst_.size() ; ++i )
            {
                step_a(i);
            }
            step_b();
            QCoreApplication::processEvents();
            #pragma omp parallel for
            for( int i=0 ; i < vvs_ptrlst_.size() ; ++i )
            {
                step_c(i);
            }
            QCoreApplication::processEvents();
            step_d();
            finish_steps();
            QCoreApplication::processEvents();
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
    void to_frame_func(const arma::mat& alpha,const arma::vec& model_func,arma::vec& frame_func);
    void to_vox_func(int index,const arma::vec& frame_func,arma::vec& vox_func);
    void proj_and_rebuild(int index,const arma::vec& vox_func,arma::vec& re_vox_func);
    void to_pix_func(int index,const arma::vec& vox_func,arma::vec& pix_func);
    void to_model_func(const arma::mat& alpha,const arma::vec& frame_func,arma::vec& model_func);
    double res_energy(const arma::vec& func0,const arma::vec& func1);
    void max_cir_cor(const arma::vec&,const arma::vec&,arma::uword& max_cir);//maximum of circular correlation
    void cir_mat(const arma::sword& cir,arma::mat& alpha);//perform cirulation on alpha
    void calc_weighted(
            const arma::fmat&vv,
            const arma::fmat&vn,
            arma::Mat<uint8_t>&vc,
            const arma::mat& alpha,
            arma::fmat&wv,
            arma::fmat&wn,
            arma::Mat<uint8_t>&wc
            );
protected:
    //iteration count
    int iter_count_;

    //max iteration number
    int max_iter_;

    //type of rotation
    RotationType rttype_;

    //input observation
    MatPtrLst vvs_ptrlst_;
    MatPtrLst vns_ptrlst_;
    CMatPtrLst vcs_ptrlst_;
    LCMatPtrLst vls_ptrlst_; //label color
    LMatPtrLst vll_ptrlst_; //labels

    //result transformation
    TsLst rt_lst_;

    //result correspondece
    DMatPtrLst alpha_ptrlst_;

    //indicating functions
    arma::mat funcs_;
    //residue function
    arma::mat res_;
    arma::vec median_res_;

    //weighted V
    MatPtrLst  wvs_ptrlst_;
    MatPtrLst  wns_ptrlst_;
    CMatPtrLst wcs_ptrlst_;

    //latent model centroid
    MatPtr xv_ptr_;
    MatPtr xn_ptr_;
    CMatPtr xc_ptr_;

    //transformed centroid
    MatPtrLst  xtv_ptrlst_;

    arma::mat var_;
    //latent model parameter
    arma::rowvec x_p_;      //probability
    arma::rowvec x_invvar_; //inversed variance

    double beta_;

    //pre-defined object info
    arma::uword obj_num_;
    arma::uvec  obj_label_;//pre-defined class label for each centroid
    std::vector<arma::uword> obj_range_;//start and end index for objects in centroid model
    arma::fmat  obj_pos_;

    //configuration
    Config::Ptr config_;
};
}
#endif // SJRCSBASE_H
