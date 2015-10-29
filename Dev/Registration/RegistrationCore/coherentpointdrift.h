#ifndef COHERENTPOINTDRIFT_H
#define COHERENTPOINTDRIFT_H
#include "registrationcore_global.h"
#include <armadillo>
#include <iostream>
#include <memory>
#include <RegistrationBase.h>
#if !defined(M_PI)
#  define M_PI 3.1415926535897932
#endif
namespace Registration {
    class REGISTRATIONCORESHARED_EXPORT CPDBase:public RegistrationBase
    {
    public:
        typedef enum{
            MaxIter,
            ErrorBelowThreshold,
            VarBelowThreshold,
            Force
        }EndMode;
        typedef struct Info{
            float omega = 0.1;// weight for uniform distribution
            int max_iter = 100;
            float fitness_th = 0.0;
            float var_th = 0.0 ;
            bool isApplyed = true;//is transform applied on input matrix source
            bool isScaled = false;
            EndMode mode;
            void* result = NULL;
        }Info;
        typedef std::shared_ptr<Info> InfoPtr;
        CPDBase():RegistrationBase(){}
        virtual ~CPDBase(){}
        //when compute in Thread
        //each col is a point
        void compute(const arma::fmat&source,const arma::fmat&target,InfoPtr&info)
        {
            reset(source,target,info);
            while(!isEnd())
            {
                computeOnce();
            }
        }

        virtual void compute()
        {
            while(!isEnd())
            {
                computeOnce();
            }
        }

        virtual void reset(const arma::fmat&source,const arma::fmat&target,InfoPtr&info)
        {
            InfoPtr_ = info;
            if(InfoPtr_->isApplyed)Y_ptr = std::shared_ptr<arma::fmat>(new arma::fmat((float*)(source.memptr()),source.n_rows,source.n_cols,false,true));
            else Y_ptr = std::shared_ptr<arma::fmat>(new arma::fmat((float*)(source.memptr()),source.n_rows,source.n_cols,true));
            X_ptr = std::shared_ptr<arma::fmat>(new arma::fmat((float*)(target.memptr()),target.n_rows,target.n_cols,false,true));
            var = 0.0;
            for( int idx = 0 ; idx < Y_ptr->n_cols ;++idx)
            {
                arma::fmat tmpa = X_ptr->each_col() - Y_ptr->col(idx);
                arma::fmat tmpb = arma::square(tmpa);
                arma::frowvec tmpc = arma::sum(tmpb);
                var += arma::mean(tmpc);
            }
            var /= Y_ptr->n_cols;
            count = 0;
        }
        void computeOnce(void)
        {
            stepE();
            stepM();
            ++count;
        }
    protected:
        virtual void stepE();
        virtual void stepM()=0;
        virtual bool isEnd()=0;
    protected:
        int count;
        float var;
        arma::fmat P_;
        std::shared_ptr<arma::fmat> X_ptr;
        std::shared_ptr<arma::fmat> Y_ptr;
        std::shared_ptr<Info> InfoPtr_;
    };

    template<typename M>
    class CPDRigid3D:public CPDBase
    {
    public:
        using CPDBase::P_;
        using CPDBase::X_ptr;
        using CPDBase::Y_ptr;
        using CPDBase::InfoPtr_;
        typedef struct Result{
            float R[9];
            float t[3];
            float s;
        }Result;
        typedef typename std::vector<typename MeshBundle<M>::Ptr> MeshList;
        CPDRigid3D();
        virtual bool initForThread(void*);
        //configure info from global configure
        bool configure(Info&);
        virtual void reset(const arma::fmat&source, const arma::fmat&target, InfoPtr &info);
        float fitness();
        virtual ~CPDRigid3D();
    protected:
        virtual void stepM();
        virtual bool isEnd();
    protected:
        float dividebyNp_;      // sum of probability
        arma::fmat A_;
        arma::fvec muX; // X center
        arma::fvec muY; // Y center
        arma::fmat mX_; // centered X
        arma::fmat mY_; // centered Y
        std::shared_ptr<Result> ResPtr_;
    };
}
#include "coherentpointdrift.hpp"
#endif // COHERENTPOINTDRIFT_H
