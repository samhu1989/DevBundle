#include "jrmpcv2.h"
namespace Registration {
template<typename M>
bool JRMPCV2<M>::configure(Config::Ptr& config,InfoPtr& info)
{
    info = std::make_shared<Info>();
    if( !config || 0==config.use_count() )
    {
        std::cerr<<"no configure"<<std::endl;
        return false;
    }else{
        if(config->has("Align_Max_Iter")){
            info->max_iter = config->getInt("Align_Max_Iter");
        }
        if(config->has("Align_Eps"))
        {
            info->eps = config->getFloat("Align_Eps");
        }
        return true;
    }
    return false;
}

template<typename M>
bool JRMPCV2<M>::initForThread(void *meshlistptr,std::vector<arma::uword>&valid_index,InfoPtr info)
{
    MeshList* list = reinterpret_cast<MeshList*>(meshlistptr);

    if(!list){
        error_string_ = "Can not locate the inputs";
        return false;
    }

    if(valid_index.size()<2){
        error_string_ = "This algorithm is designed for two multiple point clouds";
        return false;
    }

    typename MeshBundle<M>::Ptr mesh_ptr;
    std::vector<std::shared_ptr<arma::fmat>> v;
    std::vector<std::shared_ptr<arma::fmat>> vn;
    for(size_t index = 0 ; index < valid_index.size() ; ++index)
    {
        mesh_ptr = (*list)[valid_index[index]];
        float*data = (float*)mesh_ptr -> mesh_.points();
        float*ndata = (float*)mesh_ptr -> mesh_.vertex_normals();
        int N = mesh_ptr -> mesh_.n_vertices();
        if( N > 0 ){
            v.emplace_back(new arma::fmat(data,3,N,false,true));
            vn.emplace_back(new arma::fmat(ndata,3,N,false,true));
            if(!vn.back()->is_finite())
            {
                arma::uvec indices = arma::find_nonfinite(*vn.back());
                std::cerr<<"infinite in vn.back():"<<indices.size()<<"/"<<vn.back()->size()<<std::endl;
            }
        }
    }

    JRMPC<M>::initK(v,info->k);
    list->emplace_back(new MeshBundle<M>());

    while( info->k > list->back()->mesh_.n_vertices() )
    {
        list->back()->mesh_.add_vertex(typename M::Point(0,0,0));
    }

    list->back()->mesh_.request_vertex_normals();

    float* target_data = (float*)list->back()->mesh_.points();
    float* target_ndata = (float*)list->back()->mesh_.vertex_normals();
    int target_N = list->back()->mesh_.n_vertices();

    if(target_N!=info->k)
    {
        std::stringstream ss;
        ss<<"Unable to add "<<info->k<<" points to target";
        error_string_ = ss.str();
        return false;
    }

    arma::fmat x(target_data,3,target_N,false,true);
    arma::fmat xn(target_ndata,3,target_N,false,true);
    JRMPC<M>::initX(x);
    xn = arma::normalise(x);
    reset(v,vn,x,xn,info);
    JRMPC<M>::setVarColor((uint32_t*)list->back()->custom_color_.vertex_colors(),list->back()->mesh_.n_vertices());
    return true;
}

template<typename M>
void JRMPCV2<M>::reset(
        const std::vector<std::shared_ptr<arma::fmat>>&v,
        const std::vector<std::shared_ptr<arma::fmat>>&vn,
        const arma::fmat&x,
        const arma::fmat&xn,
        InfoPtr&info
        )
{
    count = 0;
    restart_count = 0;
    T_updated_ = 0;
    X_updated_ = 0;
    float* target_data = (float*)x.memptr();
    float* target_ndata = (float*)xn.memptr();
    int r = x.n_rows;
    int c = x.n_cols;
    X_ptr = std::shared_ptr<arma::fmat>(
                new arma::fmat(
                    target_data,
                    r,
                    c,
                    false,
                    true
                    )
            );
    Xn_ptrs = std::shared_ptr<arma::fmat>(
                new arma::fmat(
                    target_ndata,
                    r,
                    c,
                    false,
                    true
                    )
            );
    V_ptrs = v;
    Vn_ptrs = vn;
    info_ptr = info;
    P_ = arma::frowvec(c);
    P_.fill(1.0/float(c));
    //initialize R and t
    res_ptr = ResPtr(new Result);
    std::vector<std::shared_ptr<arma::fmat>>::iterator VIter;
    std::vector<std::shared_ptr<arma::fmat>>::iterator VnIter;
    arma::fvec meanX = arma::mean(*X_ptr,1);
    arma::fmat& v0 = *V_ptrs.front();
    arma::fmat& vn0 = *Vn_ptrs.front();
    VnIter = Vn_ptrs.begin();
    for( VIter = V_ptrs.begin() ; VIter != V_ptrs.end() ; ++VIter )
    {
        arma::fmat&v = **VIter;
        arma::fmat&vn = **VnIter;

//        arma::fvec rand(v.n_cols/2+1,arma::fill::randu);
//        arma::uvec select = arma::conv_to<arma::uvec>::from(v.n_cols*rand);
//        arma::fvec meanV = arma::mean(v.cols(select),1);
        arma::fvec meanV = arma::mean(v,1);
        res_ptr->Rs.emplace_back(new arma::fmat(3,3,arma::fill::eye));
        res_ptr->ts.emplace_back(new arma::fvec( meanX - meanV ));
        v.each_col() += *(res_ptr->ts.back());
        if(VIter!=V_ptrs.begin())
        {
            alignAroundZ(v0,vn0,v,vn,*res_ptr->Rs.back());
            v = (*res_ptr->Rs.back())*v;
            vn = (*res_ptr->Rs.back())*vn;
            vn = arma::normalise(vn);
            (*res_ptr->ts.back()) = *res_ptr->Rs.back()*(*res_ptr->ts.back());
        }
        ++VnIter;
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

    var = arma::frowvec(X_ptr->n_cols);
    var.fill(1.0/maxvar);

    X_sum = arma::fmat(X_ptr->n_rows,X_ptr->n_cols,arma::fill::zeros);
    var_sum = arma::frowvec(X_ptr->n_cols,arma::fill::zeros);
    alpha_sum = arma::frowvec(X_ptr->n_cols,arma::fill::zeros);
    alpha_sumij = arma::frowvec(X_ptr->n_cols,arma::fill::zeros);

    float h;
    h = 2.0 / arma::mean(var);
    float gamma = info_ptr->gamma;
    beta = gamma/(h*(gamma+1));

    std::cerr<<V_ptrs.size()<<std::endl;
    std::cerr<<"V_pts[0]:"<<V_ptrs[0]->n_cols<<std::endl;
    std::cerr<<"V_pts[1]:"<<V_ptrs[1]->n_cols<<std::endl;
    std::cerr<<"X_ptr->n_cols:"<<X_ptr->n_cols<<std::endl;
}


}

