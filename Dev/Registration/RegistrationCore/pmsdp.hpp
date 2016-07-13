#ifndef PMSDP_HPP
#define PMSDP_HPP
#include "pmsdp.h"
namespace Registration {
template<typename M>
PMSDP<M>::PMSDP()
{

}
template<typename M>
PMSDP<M>::~PMSDP()
{

}

template<typename M>
bool PMSDP<M>::configure(Config::Ptr&,InfoPtr&)
{
    std::cerr<<"PMSDP::configure(Config::Ptr&,InfoPtr&)"<<std::endl;
    std::cerr<<"Under Implementation"<<std::endl;
    return false;
}

template<typename M>
void PMSDP<M>::generateObj()
{
    std::cerr<<"PMSDP::generateObj()"<<std::endl;
    std::cerr<<"Under Implementation"<<std::endl;
}

template<typename M>
void PMSDP<M>::generateConstraint()
{
    std::cerr<<"PMSDP::generateConstraint()"<<std::endl;
    std::cerr<<"Under Implementation"<<std::endl;
}
template<typename M>
void PMSDP<M>::projectRX()
{
    std::cerr<<"PMSDP::projectRX()"<<std::endl;
    std::cerr<<"Under Implementation"<<std::endl;
}
}
#endif // PMSDP_HPP
