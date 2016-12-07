#ifndef SJRCSBASE_H
#define SJRCSBASE_H
#include "jrcscore_global.h"
#include "common.h"
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
    typedef struct{
        float R[9];
        float t[3];
    }T;
    typedef std::vector<T> Ts;
    typedef std::vector<Ts> TsLst;
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
    virtual void compute(void)
    {
        while(!isEnd())
        {
            prepare_compute();
            #pragma omp parallel for
            for( int i=0 ; i < vvs_ptrlst_.size() ; ++i )
            {
                parallel_compute(i);
            }
            finish_compute();
        }
    }
protected:
    virtual bool isEnd(void);
    virtual void prepare_compute(void);
    virtual void parallel_compute(int i);
    virtual void finish_compute(void);
protected：
	//Calculate Alpha(Expectation)
	//Construct Residue Function
	
	////Median of Residue Function
	
	//Calculate Correlation Function
	//Alter Alpha
		
	//Update RT
	
	////Update X
	
	////Update var pk
protected:
    //iteration count
    int iter_count_;
    //max iteration number
    int max_iter_;
    //input observation
    MatPtrLst vvs_ptrlst_;
    MatPtrLst vns_ptrlst_;
    CMatPtrLst vcs_ptrlst_;
    LCMatPtrLst vls_ptrlst_; //label color
    LMatPtrLst vll_ptrlst_; //labels
    //result transformation
    TsLst rt_lst_;

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

    //sub-view point to same memory space of latent model centroids
    MatPtrLst   objv_ptrlst_;
    MatPtrLst   objn_ptrlst_;
    CMatPtrLst  objc_ptrlst_;

    //configuration
    Config::Ptr config_;
};
}
#endif // SJRCSBASE_H
