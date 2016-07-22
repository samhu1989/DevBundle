#ifndef PMSDP_H
#define PMSDP_H
#include "RegistrationBase.h"
#include <armadillo>
#include "common.h"
#include "sdp.h"
#include "pmsdp_matlab.h"
namespace Registration{
template<typename M>
class PMSDP:public RegistrationBase{
    using RegistrationBase::end_;
public:
    typedef typename std::vector<typename MeshBundle<M>::Ptr> MeshList;
    typedef struct Result{
        float R[9];
        float t[3];
    }Result;
    typedef std::shared_ptr<Result> ResPtr;
    typedef struct Info{
        void* result = NULL;
    }Info;
    typedef std::shared_ptr<Info> InfoPtr;
    PMSDP();
    virtual ~PMSDP();
    bool configure(Config::Ptr&,InfoPtr&);
    virtual bool initForThread(void*,InfoPtr&);
    virtual bool initForThread(void*,std::vector<arma::uword>&,InfoPtr&){
        std::cerr<<"PMSDP::initForThread(void*,std::vector<arma::uword>&,InfoPtr&)"<<std::endl;
        std::cerr<<"Under Implementation"<<std::endl;
        return false;
    }
    virtual void compute(void)
    {
        generateObj();
        generateConstraint();
        if(!sdp_.init())
        {
            computeByMatlab();
        }else{
            assert(sdp_.solve());
            projectRX();
        }
    }
    //configure info from global configure
    ResPtr result(void){return res_ptr_;}
protected:
    virtual void computeByMatlab(void);
    virtual void generateObj(void);
    virtual void generateConstraint(void);
    virtual void projectRX(void);
private:
    Optimization::SDP sdp_;
    ResPtr res_ptr_;
    std::shared_ptr<arma::fmat> P_;
    std::shared_ptr<arma::fmat> Q_;
    std::shared_ptr<PMSDP_MATLAB> pmsdp_matlab_;
};
}
#endif // PMSDP_H
