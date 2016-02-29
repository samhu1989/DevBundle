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
    using typename JRMPC<M>::ResPtr;
    JRMPCV2():JRMPC<M>(){;}
    virtual bool configure(Config::Ptr&,InfoPtr&);
    virtual ResPtr result(void)
    {
        std::cerr<<"JRMPCV2 Result"<<std::endl;
        return res_ptr;
    }
    virtual void compute(void)
    {
        std::cerr<<"JRMPCV2 compute"<<std::endl;
        do{
            JRMPC<M>::computeOnce();
            JRMPC<M>::varToColor();
            ++count;
        }while(!JRMPC<M>::isEnd());
        res_ptr->X = X_ptr;
    }
protected:
    std::vector<std::shared_ptr<arma::fmat>> VN_ptrs;
    using JRMPC<M>::res_ptr;
    using JRMPC<M>::X_ptr;
    using JRMPC<M>::count;
};

}
#include "jrmpcv2.hpp"
#endif // LJRMPC_H
