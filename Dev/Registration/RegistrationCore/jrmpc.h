#ifndef JRMPC_H
#define JRMPC_H
#include "RegistrationBase.h"
#include <armadillo>
namespace Registration {
    template<typename M>
    class JRMPC:public RegistrationBase
    {
        typedef struct Info{

        }Info;

        typedef typename std::vector<typename MeshBundle<M>::Ptr> MeshList;
        typedef std::shared_ptr<Info> InfoPtr;
    public:
        virtual void initForThread(void *meshlistptr);

        virtual void reset(const arma::fmat&source,const arma::fmat&target,InfoPtr&info);

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
    private:
        MeshList mesh_list;
    };
}
#include <jrmpc.hpp>
#endif // JRMPC_H
