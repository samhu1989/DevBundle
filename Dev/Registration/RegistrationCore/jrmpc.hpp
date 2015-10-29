#ifndef JRMPC_HPP
#define JRMPC_HPP
#include <memory>
#include <strstream>
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

    std::shared_ptr<Info> info(new Info);
    if(!configure(*info))
    {
        error_string_ = "With no proper configure";
        return false;
    }

    typename MeshList::iterator meshIter;
    std::vector<std::shared_ptr<arma::fmat>> source;
    for(meshIter=list->begin();meshIter!=list->end();++meshIter)
    {
        float*data = (float*)(*meshIter) -> mesh_.points();
        int N = (*meshIter) -> mesh_.n_vertices();
        source.push_back(
                    std::shared_ptr<arma::fmat>(new arma::fmat(data,3,N,false,true))
                    );
    }

    initK(source,info->k);
    list->push_back(typename MeshBundle<M>::Ptr(new MeshBundle<M>()));

    while( info->k > list->back()->mesh_.n_vertices() )
    {
        list->back()->mesh_.add_vertex(typename M::Point(0,0,0));
    }

    float* target_data = (float*)list->back()->mesh_.points();
    int target_N = list->back()->mesh_.n_vertices();

    if(target_N!=info->k)
    {
        std::stringstream ss;
        ss<<"Unable to add "<<info->k<<" points to target";
        error_string_ = ss.str();
        return false;
    }

    arma::fmat target(target_data,3,target_N,false,true);
    initX(source,target);
    reset(source,target,info);
    return true;
}

template<typename M>
void JRMPC<M>::reset(
        const std::vector<std::shared_ptr<arma::fmat>>&source,
        const arma::fmat&target,
        InfoPtr&info
        )
{
    count = 0;
    float* target_data = (float*)target.memptr();
    int r = target.n_rows;
    int c = target.n_cols;
    X_ptr = std::shared_ptr<arma::fmat>(
                new arma::fmat(
                    target_data,
                    r,
                    c,
                    false,
                    true
                    )
            );
    V_ptrs = source;
    info_ptr = info;
    P_ = arma::fvec(c);
    P_.fill(1.0/float(c));
    //initialize R and t
    res_ptr = ResPtr(new Result);
    std::vector<std::shared_ptr<arma::fmat>>::iterator VIter;
    arma::fvec meanX = arma::mean(*X_ptr,1);
    for( VIter = V_ptrs.begin() ; VIter != V_ptrs.end() ; ++VIter )
    {
        arma::fmat&v = **VIter;
        arma::fvec meanV = arma::mean(v,1);
        res_ptr->Rs.push_back(std::make_shared<arma::fmat>(3,3,arma::fill::eye));
        res_ptr->ts.push_back(std::make_shared<arma::fvec>(meanX-meanV));
        v.each_col() += *(res_ptr->ts.back());
    }
    info_ptr ->result = (void*)(res_ptr.get());

    arma::fvec maxAllXYZ,minAllXYZ;

    for( VIter = V_ptrs.begin() ; VIter != V_ptrs.end() ; ++VIter )
    {
        arma::fmat&v = **VIter;
        arma::fvec maxVXYZ = arma::max(v,1);
        arma::fvec minVXYZ = arma::min(v,1);
        if(VIter==V_ptrs.begin())
        {
            maxAllXYZ = maxVXYZ;
            minAllXYZ = minVXYZ;
        }else{
            maxAllXYZ = arma::max(maxVXYZ,maxAllXYZ);
            minAllXYZ = arma::min(minVXYZ,minAllXYZ);
        }
    }

    float maxvar = arma::accu(arma::square(maxAllXYZ-minAllXYZ));

    (*X_ptr)*=(0.5*std::sqrt(maxvar));

    var = arma::fvec(X_ptr->n_cols);
    var.fill(1.0/maxvar);

    float h;
    h = 2 / arma::mean(var);
    float gamma = info_ptr->gamma;
    beta = gamma/(h*(gamma+1));
}

template<typename M>
void JRMPC<M>::initK(
        const std::vector<std::shared_ptr<arma::fmat>>&source,
        int&k
        )
{
    //if the k is already set during configure
    if(0!=k)return;
    std::vector<std::shared_ptr<arma::fmat>>::const_iterator iter;
    arma::fvec nviews(source.size());
    int idx = 0;
    for( iter = source.begin() ; iter != source.end()  ; ++iter )
    {
        nviews(idx) = (*iter)->n_cols;
        ++idx;
    }
    k = int(0.8*arma::median(nviews));
    k = k > 12 ? k : 12 ;
}

template<typename M>
void JRMPC<M>::initX(const std::vector<std::shared_ptr<arma::fmat>>&source,
        arma::fmat&target
        )
{
    int k = target.n_cols;
    arma::frowvec az = arma::randu<arma::frowvec>(k);
    arma::frowvec el = arma::randu<arma::frowvec>(k);
    az*=2*M_PI;
    el*=2*M_PI;
    target.row(0) = arma::cos(az)%arma::cos(el);
    target.row(1) = arma::sin(el);
    target.row(2) = arma::sin(az)%arma::cos(el);

}

template<typename M>
void JRMPC<M>::stepE()
{
    ;
}

template<typename M>
void JRMPC<M>::stepMa()
{
    ;
}

template<typename M>
void JRMPC<M>::stepMb()
{
    ;
}

template<typename M>
void JRMPC<M>::stepMc()
{
    ;
}

template<typename M>
void JRMPC<M>::stepMd()
{
    ;
}

template<typename M>
bool JRMPC<M>::isEnd()
{
    if( count >= info_ptr->max_iter )
    {
        info_ptr->mode = MaxIter;
        return true;
    }
    return false;
}

}
#endif

