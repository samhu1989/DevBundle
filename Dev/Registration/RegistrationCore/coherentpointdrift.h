#ifndef COHERENTPOINTDRIFT_H
#define COHERENTPOINTDRIFT_H
#include <armadillo>
#include <iostream>
#include <memory>
#if !defined(M_PI)
#  define M_PI 3.1415926535897932
#endif
namespace Registration {
    class CPDBase
    {
    public:
        typedef enum{
            MaxIter,
            ErrorBelowThreshold,
            VarBelowThreshold
        }EndMode;
        typedef struct Info{
            float omega = 0.3;// weight for uniform distribution
            int max_iter = 10;
            float fitness_th = std::numeric_limits<float>::max();
            float var_th = std::numeric_limits<float>::max() ;
            bool isApplyed = true;//is transform applied on input matrix source
            bool isScaled = false;
            EndMode mode;
            void* result = NULL;
        }Info;
        CPDBase(){}
        virtual ~CPDBase(){}
        //each col is a point
        void compute(const arma::fmat&source,const arma::fmat&target,Info&info)
        {
            reset(source,target,info);
            while(!isEnd())
            {
                computeOnce();
            }
        }
        virtual void reset(const arma::fmat&source,const arma::fmat&target,Info&info)
        {
            InfoPtr_ = std::shared_ptr<Info>(&info);
            if(InfoPtr_->isApplyed)Y_ = arma::fmat((float*)(source.memptr()),source.n_rows,source.n_cols,false,true);
            else Y_ = arma::fmat((float*)(source.memptr()),source.n_rows,source.n_cols,true);
            X_ = arma::fmat((float*)(target.memptr()),target.n_rows,target.n_cols,false,true);
            var = 0.0;
            for( int idx = 0 ; idx < Y_.n_cols ;++idx)
            {
                arma::fmat tmpm = X_.each_col() - Y_.col(idx);
                arma::fvec tmpv = arma::sum(arma::square(tmpm));
                var += arma::mean(tmpv);
            }
            var /= Y_.n_cols;
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
        arma::fmat X_;
        arma::fmat Y_;
        std::shared_ptr<Info> InfoPtr_;
    };

    class CPDRigid3D:public CPDBase
    {
    public:
        using CPDBase::P_;
        using CPDBase::X_;
        using CPDBase::Y_;
        using CPDBase::InfoPtr_;
        typedef struct Result{
            float R[9];
            float t[3];
            float s;
        }Result;
        CPDRigid3D();
        virtual void reset(const arma::fmat&source,const arma::fmat&target,Info&info);
        float fitness();
        virtual ~CPDRigid3D();
    protected:
        virtual void stepM();
        virtual bool isEnd();
        virtual void updateCenter();
        virtual void updateRotation();
        virtual void updateScale();
        virtual void updateTranslation();
        virtual void updateVar();
        virtual void updateY();
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

#endif // COHERENTPOINTDRIFT_H
