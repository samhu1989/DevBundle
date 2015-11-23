#ifndef REGISTRATIONTHREAD_HPP
#define REGISTRATIONTHREAD_HPP
#include "RegistrationThreadT.h"
namespace Registration{

template<typename Reg,typename M>
RegistrationThreadT<Reg,M>::RegistrationThreadT(QObject* parent):QThread(parent)
{
    setTerminationEnabled(true);
}

template<typename Reg,typename M>
bool RegistrationThreadT<Reg,M>::init(MeshList& mesh_list,Config::Ptr& config)
{
    typename Reg::InfoPtr info;
    reg_.configure(config,info);
    return reg_.initForThread((void*)&mesh_list,info);
}

template<typename Reg,typename M>
bool RegistrationThreadT<Reg,M>::init(MeshList& mesh_list,std::vector<arma::uword>& valid_index,Config::Ptr& config)
{
    typename Reg::InfoPtr info;
    if(!reg_.configure(config,info))return false;
    return reg_.initForThread((void*)&mesh_list,valid_index,info);
}

template<typename Reg,typename M>
void RegistrationThreadT<Reg,M>::compute(void)
{
    reg_.compute();
}

template<typename Reg,typename M>
typename Reg::ResPtr RegistrationThreadT<Reg,M>::result(void)
{
    return reg_.result();
}

}
#endif
