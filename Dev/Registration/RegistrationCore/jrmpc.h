#ifndef JRMPC_H
#define JRMPC_H
#include "RegistrationBase.h"
#include <armadillo>
namespace Registration {
    template<typename M>
    class JRMPC:public RegistrationBase
    {
        typedef enum{
            MaxIter,
            ErrorBelowThreshold,
            VarBelowThreshold
        }EndMode;

        typedef struct Info{
            float gama = 0.2;// weight for uniform distribution
            int max_iter = 50;
            float fitness_th = 0.0;
            float var_th = 0.0 ;
            bool isApplyed = true;//is transform applied on input matrix source
            EndMode mode;
            void* result = NULL;
        }Info;

        typedef struct Result{
            std::shared_ptr<arma::fmat> Rs;
            std::shared_ptr<arma::fvec> ts;
        }Result;

        typedef typename std::vector<typename MeshBundle<M>::Ptr> MeshList;
        typedef std::shared_ptr<Info> InfoPtr;
    public:
        JRMPC();
        bool configure(Info&);
        virtual bool initForThread(void *meshlistptr);

        virtual void reset(
                const std::vector<std::shared_ptr<arma::fmat>>&source,
                const arma::fmat&target,
                InfoPtr&info
                );

        virtual void compute(void)
        {
            while(!isEnd())
            {
                computeOnce();
            }
        }

        virtual void computeOnce(void)
        {
            stepE();
            stepMa();
            stepMb();
            stepMc();
            stepMd();
            ++count;
        }

    protected:
        virtual void stepE();
        virtual void stepMa();
        virtual void stepMb();
        virtual void stepMc();
        virtual void stepMd();
        virtual bool isEnd();

    protected:
        int count;
        arma::fmat P_;

        std::shared_ptr<arma::fmat> X_ptr;
        arma::fvec var;

        std::vector<std::shared_ptr<arma::fmat>> V_ptrs;
        std::vector<std::shared_ptr<arma::fmat>> alpha_ptrs;
        InfoPtr info_ptr;
    };
}
#include <jrmpc.hpp>
#endif // JRMPC_H
