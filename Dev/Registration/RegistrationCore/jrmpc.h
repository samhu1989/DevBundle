#ifndef JRMPC_H
#define JRMPC_H
#include "RegistrationBase.h"
#include <armadillo>
#include "common.h"
namespace Registration {
    template<typename M>
    class JRMPC:public RegistrationBase
    {
        using RegistrationBase::end_;
    public:
        typedef enum{
            MaxIter,
            Converge,
            ErrorBelowThreshold,
            VarBelowThreshold,
            Force
        }EndMode;

        typedef struct Info{
            int k = 0;
            float gamma = 0.05;// weight for uniform distribution
            int max_iter = 50;
            float fitness_th = 0.0;
            float var_th = 0.0 ;
            float eps = 1e-7;
            bool isApplyed = true;//is transform applied on input matrix source
            EndMode mode;
            void* result = NULL;
        }Info;

        typedef struct Result{
            std::vector<std::shared_ptr<arma::fmat>> Rs;
            std::vector<std::shared_ptr<arma::fvec>> ts;
            std::shared_ptr<arma::fmat> X;
        }Result;

        typedef typename std::vector<typename MeshBundle<M>::Ptr> MeshList;
        typedef std::shared_ptr<Info> InfoPtr;
        typedef std::shared_ptr<Result> ResPtr;
        typedef std::shared_ptr<arma::fmat> MatPtr;
    public:
        JRMPC();
        virtual bool configure(Config::Ptr&,InfoPtr&);
        virtual ResPtr result(void){return res_ptr;}

        virtual bool initForThread(void *meshlistptr,InfoPtr info);
        virtual bool initForThread(void *meshlistptr,std::vector<arma::uword>&valid_index,InfoPtr info);

        virtual void reset(
                const std::vector<std::shared_ptr<arma::fmat> > &source,
                InfoPtr &info
                )
        {
            arma::fmat target;
            initK(source,info->k);
            initX(source,target);
            reset(source,target,info);
        }

        virtual void reset(
                const std::vector<std::shared_ptr<arma::fmat>>&source,
                const arma::fmat&target,
                InfoPtr&info
                );

        virtual void compute(
                const std::vector<std::shared_ptr<arma::fmat>>&source,
                InfoPtr&info)
        {
            reset(source,info);
            compute();
        }

        virtual void compute(
                const std::vector<std::shared_ptr<arma::fmat>>&source,
                const arma::fmat&target,
                InfoPtr&info)
        {
            reset(source,target,info);
            compute();
        }

        virtual void compute(void)
        {

            do{
                computeOnce();
                varToColor();
                ++count;
            }while(!isEnd());
            res_ptr->X = X_ptr;
        }

        virtual void computeOnceStep(void)
        {
            stepE();
            stepMa();
            stepMbc();
            stepMd();
        }

        virtual void computeOnce(void);
        virtual void setVarColor(uint32_t*,int k);

    protected:
        virtual void initK(const std::vector<std::shared_ptr<arma::fmat>>&source,int&k);
        virtual void initX(arma::fmat&target);
        virtual void initX(const std::vector<std::shared_ptr<arma::fmat>>&source,arma::fmat&target);
        virtual void stepE();
        virtual void stepMa();
        virtual void stepMbc();
        virtual void stepMd();
        virtual bool isEnd();
        virtual void varToColor();
    protected:
        int count;

        arma::frowvec P_;
        std::shared_ptr<arma::fmat> X_ptr;
        arma::frowvec var;
        std::shared_ptr<arma::Col<uint32_t>> var_color;

        std::vector<std::shared_ptr<arma::fmat>> V_ptrs;
        std::vector<std::shared_ptr<arma::fmat>> alpha_ptrs;
        InfoPtr info_ptr;
        ResPtr res_ptr;

        float beta;

        arma::fmat X_sum;
        arma::frowvec var_sum;
        arma::frowvec alpha_sum;
        arma::frowvec alpha_sumij;

        int T_updated_;
        int X_updated_;
    };
}
#include <jrmpc.hpp>
#endif // JRMPC_H
