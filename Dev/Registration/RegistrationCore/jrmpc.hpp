#ifndef JRMPC_HPP
#define JRMPC_HPP
#include "jrmpc.h"
namespace Registration {
template<typename M>
JRMPC<M>::JRMPC():
    count(0)
{

}

template<typename M>
bool JRMPC<M>::configure(Info&)
{
    return true;
}

template<typename M>
bool JRMPC<M>::initForThread(void *meshlistptr)
{
    MeshList* list = reinterpret_cast<MeshList*>(meshlistptr);
    if(!list){
        error_string_ = "Can not locate the inputs";
        return false;
    }
    if(list->size()<2){
        error_string_ = "This algorithm is designed for two multiple point clouds";
        return false;
    }
    arma::fmat source((float*)(*list)[0]->mesh_.points(),3,(*list)[0]->mesh_.n_vertices(),false,true);
    arma::fmat target((float*)(*list)[1]->mesh_.points(),3,(*list)[1]->mesh_.n_vertices(),false,true);

    std::shared_ptr<Info> info(new Info);
    if(!configure(*info))
    {
        error_string_ = "With no proper configure";
        return false;
    }

//  reset(source,target,info);
    return true;
}

template<typename M>
void JRMPC<M>::reset(
        const std::vector<std::shared_ptr<arma::fmat>>&source,
        const arma::fmat&target,
        InfoPtr&info
        )
{
    ;
}

template<typename M>
void JRMPC<M>::stepE()
{
    ;
}

template<typename M>
void JRMPC<M>::stepMa()
{

}

template<typename M>
void JRMPC<M>::stepMb()
{

}

template<typename M>
void JRMPC<M>::stepMc()
{

}

template<typename M>
void JRMPC<M>::stepMd()
{

}

template<typename M>
bool JRMPC<M>::isEnd()
{
    ;
}

}
#endif

