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
bool RegistrationThreadT<Reg,M>::init(MeshList& mesh_list)
{
    return reg_.initForThread((void*)&mesh_list);
}

template<typename Reg,typename M>
void RegistrationThreadT<Reg,M>::compute(void)
{
    reg_.compute();
}

}
#endif
