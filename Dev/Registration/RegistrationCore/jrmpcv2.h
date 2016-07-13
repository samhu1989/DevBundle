#ifndef LJRMPC_H
#define LJRMPC_H
#include "jrmpc.h"
namespace Registration {

template<typename M>
class JRMPCV2:public JRMPC<M>
{
public:
    using typename JRMPC<M>::Info;
    using typename JRMPC<M>::InfoPtr;
    using typename JRMPC<M>::Result;
    using typename JRMPC<M>::ResPtr;
    using typename JRMPC<M>::MeshList;
    JRMPCV2():JRMPC<M>(){;}
    virtual bool configure(Config::Ptr&,InfoPtr&);
    virtual bool initForThread(
            void *meshlistptr,
            std::vector<arma::uword>&valid_index,
            InfoPtr info
            );
    virtual void reset(
            const std::vector<std::shared_ptr<arma::fmat>>&v,
            const std::vector<std::shared_ptr<arma::fmat>>&vn,
            const arma::fmat&x,
            const arma::fmat&xn,
            InfoPtr&info
            );
protected:
    std::vector<std::shared_ptr<arma::fmat>> Vn_ptrs;
    std::shared_ptr<arma::fmat> Xn_ptrs;
    using JRMPC<M>::alignAroundZ;
    using JRMPC<M>::info_ptr;
    using JRMPC<M>::res_ptr;
    using JRMPC<M>::X_ptr;
    using JRMPC<M>::count;
    using JRMPC<M>::restart_count;
    using JRMPC<M>::T_updated_;
    using JRMPC<M>::X_updated_;
    using JRMPC<M>::V_ptrs;
    using JRMPC<M>::P_;
    using JRMPC<M>::X_sum;
    using JRMPC<M>::var;
    using JRMPC<M>::var_sum;
    using JRMPC<M>::alpha_sum;
    using JRMPC<M>::alpha_sumij;
    using JRMPC<M>::beta;
    using JRMPC<M>::error_string_;
};

}
#endif // LJRMPC_H
