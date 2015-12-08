#ifndef LJRMPC_H
#define LJRMPC_H
#include "jrmpc.h"
namespace Registration {

template<typename M>
class LJRMPC:public JRMPC<M>
{
public:
    LJRMPC():JRMPC<M>(){;}
    virtual bool configure(Config::Ptr&,InfoPtr&);
    virtual ResPtr result(void){return res_ptr;}
signals:

private:
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
#include "ljrmpc.hpp"
#endif // LJRMPC_H
