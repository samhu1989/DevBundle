#ifndef JRMPC_HPP
#define JRMPC_HPP
#include <memory>
#include <strstream>
#include "jrmpc.h"
#include "common.h"
#include <iomanip>
namespace Registration {
template<typename M>
JRMPC<M>::JRMPC():RegistrationBase(),
    count(0)
{

}

template<typename M>
bool JRMPC<M>::configure(Config::Ptr& config,InfoPtr& info)
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
            info->eps_ = config->getFloat("Align_Eps");
        }
        return true;
    }
    return false;
}

template<typename M>
bool JRMPC<M>::initForThread(void *meshlistptr,InfoPtr info)
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

    typename MeshList::iterator meshIter;
    std::vector<std::shared_ptr<arma::fmat>> source;
    for(meshIter=list->begin();meshIter!=list->end();++meshIter)
    {
        float*data = (float*)(*meshIter) -> mesh_.points();
        int N = (*meshIter) -> mesh_.n_vertices();
        source.emplace_back(new arma::fmat(data,3,N,false,true));
    }

    initK(source,info->k);
    list->emplace_back(new MeshBundle<M>());

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
    initX(target);
    reset(source,target,info);
    setVarColor((uint32_t*)list->back()->custom_color_.vertex_colors(),list->back()->mesh_.n_vertices());
    return true;
}

template<typename M>
bool JRMPC<M>::initForThread(void *meshlistptr,std::vector<arma::uword>&valid_index,InfoPtr info)
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
    std::vector<std::shared_ptr<arma::fmat>> source;
    for(size_t index = 0 ; index < valid_index.size() ; ++index)
    {
        mesh_ptr = (*list)[valid_index[index]];
        float*data = (float*)mesh_ptr -> mesh_.points();
        int N = mesh_ptr -> mesh_.n_vertices();
        if( N > 0 )source.emplace_back(new arma::fmat(data,3,N,false,true));
    }

    initK(source,info->k);
    list->emplace_back(new MeshBundle<M>());

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
    initX(target);
    reset(source,target,info);
    setVarColor((uint32_t*)list->back()->custom_color_.vertex_colors(),list->back()->mesh_.n_vertices());
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
    restart_count = 0;
    T_updated_ = 0;
    X_updated_ = 0;
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
    P_ = arma::frowvec(c);
    P_.fill(1.0/float(c+1.0));
    //initialize R and t
    res_ptr = ResPtr(new Result);
    std::vector<std::shared_ptr<arma::fmat>>::iterator VIter;
    arma::fvec meanX = arma::mean(*X_ptr,1);
    for( VIter = V_ptrs.begin() ; VIter != V_ptrs.end() ; ++VIter )
    {
        arma::fmat&v = **VIter;
//        arma::fvec rand(v.n_cols/2+1,arma::fill::randu);
//        arma::uvec select = arma::conv_to<arma::uvec>::from(v.n_cols*rand);
//        arma::fvec meanV = arma::mean(v.cols(select),1);
        arma::fvec meanV = arma::mean(v,1);
        res_ptr->Rs.emplace_back(new arma::fmat(3,3,arma::fill::eye));
        res_ptr->ts.emplace_back(new arma::fvec( meanX - meanV ));
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
    std::cerr<<"beta:"<<beta<<std::endl;
    std::cerr<<"V_pts[0]:"<<V_ptrs[0]->n_cols<<std::endl;
    std::cerr<<"V_pts[1]:"<<V_ptrs[1]->n_cols<<std::endl;
    std::cerr<<"X_ptr->n_cols:"<<X_ptr->n_cols<<std::endl;
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
    k = int(0.5*arma::median(nviews));
    k = k > 300 ? k : 300 ;
    k = k < 8000 ? k : 8000;
}

template<typename M>
void JRMPC<M>::initX(arma::fmat&target)
{
    int k = target.n_cols;
    arma::frowvec az = arma::randu<arma::frowvec>(k);
    arma::frowvec el = arma::randu<arma::frowvec>(k);
    az*=2*M_PI;
    el*=2*M_PI;
    target.row(0) = 0.4*arma::cos(az)%arma::cos(el);
    target.row(1) = 0.4*arma::sin(el);
    target.row(2) = 0.4*arma::sin(az)%arma::cos(el);
}

template<typename M>
void JRMPC<M>::initX(const std::vector<std::shared_ptr<arma::fmat>>&source,
        arma::fmat&target
        )
{
    typedef std::vector<std::shared_ptr<arma::fmat>> MatPtrList;
    int k = target.n_cols;
    target.fill(0.0);
    arma::uvec rand_indice;
    arma::vec rand;
    MatPtrList::const_iterator iter;
    float cnt = 0.0;
    for(iter=source.cbegin();iter!=source.cend();++iter)
    {
        if( (*iter)->n_cols > k )
        {
            rand = arma::randu<arma::vec>(k);
            rand *= double( (*iter)->n_cols );
            rand_indice = arma::conv_to<arma::uvec>::from(rand);
            target += 0.95*(**iter).cols(rand_indice);
            target += 0.05*arma::randn<arma::fmat>(3,k);
            cnt += 1.0;
        }
    }
    target /= cnt;
}

template<typename M>
void JRMPC<M>::stepE()
{
    arma::fmat& X_ = *X_ptr;
    std::vector<MatPtr>::iterator VIter;
    int idx = 0;
    for(VIter=V_ptrs.begin();VIter!=V_ptrs.end();++VIter)
    {
        arma::fmat& V_ = **VIter;
        while(alpha_ptrs.size()<=idx)
        {
            alpha_ptrs.emplace_back(new arma::fmat(V_.n_cols,X_.n_cols));
        }
        arma::fmat& alpha = *alpha_ptrs[idx];
        for(int r=0;r<alpha.n_rows;++r)
        {
            arma::fmat tmpm = X_.each_col() - V_.col(r);
            alpha.row(r) = arma::sum(arma::square(tmpm));
        }
        alpha.each_row()%=(-0.5*var);
        alpha = arma::exp(alpha);
        alpha.each_row()%=arma::pow(var,1.5);
        alpha.each_row()%=P_.t();
        arma::fvec alpha_rowsum = arma::sum(alpha,1)+beta;
        alpha.each_col()/=alpha_rowsum;
        ++idx;
    }
}

template<typename M>
void JRMPC<M>::stepMa()
{
    arma::fmat& X_ = *X_ptr;
    std::vector<MatPtr>::iterator VIter;
    int idx = 0;
    for(VIter=V_ptrs.begin();VIter!=V_ptrs.end();++VIter)
    {
        arma::fmat& V_ = **VIter;
        arma::fmat& R = *(res_ptr->Rs[idx]);
        arma::fmat dR;
        arma::fvec& t = *(res_ptr->ts[idx]);
        arma::fvec dt;
        arma::fmat& alpha = *alpha_ptrs[idx];
        arma::frowvec alpha_colsum = arma::sum(alpha);
        arma::frowvec square_lambda = (var)%alpha_colsum;
        arma::fmat W = V_*alpha;
        W.each_row()/=alpha_colsum;
        arma::frowvec p(square_lambda.n_cols,arma::fill::ones);
        arma::frowvec square_norm_lambda = square_lambda / arma::accu(square_lambda);
        p -= square_norm_lambda;
        arma::fmat tmp = X_;
        tmp.each_row() %= p;
        tmp.each_row() %= square_lambda;
        arma::fmat A = tmp*(W.t());
        arma::fmat U,V;
        arma::fvec s;
        if(!arma::svd(U,s,V,A,"std"))
        {
            A += 1e-6*arma::fmat(3,3,arma::fill::eye);
            if(!arma::svd(U,s,V,A))
            {
                U = arma::fmat(3,3,arma::fill::eye);
                V = U;
                s = arma::fvec(3,arma::fill::ones);
                end_ = true;
            }
        }
        arma::fmat C(3,3,arma::fill::eye);
        C(2,2) = arma::det( U * V.t() )>=0 ? 1.0 : -1.0;
        dR = U*C*(V.t());
        R = dR*R;
        tmp = X_ - dR*W;
//        tmp.each_row()%=square_norm_lambda;
//        dt = arma::sum(tmp,1);
        dt = arma::mean(tmp,1);
        t = dR*t + dt;
        V_ = dR*V_;
        V_.each_col() += dt;
        ++idx;
    }
}

template<typename M>
void JRMPC<M>::stepMbc()
{
    arma::fmat& X_ = *X_ptr;
    arma::fmat X_sum(X_.n_rows,X_.n_cols,arma::fill::zeros);
    arma::frowvec var_sum(X_.n_cols,arma::fill::zeros);
    arma::frowvec alpha_sum(X_.n_cols,arma::fill::zeros);
    std::vector<MatPtr>::iterator VIter;
    int idx = 0;
    for(VIter=V_ptrs.begin();VIter!=V_ptrs.end();++VIter)
    {
        arma::fmat& V_ = **VIter;
        arma::fmat& alpha = *alpha_ptrs[idx];
        arma::frowvec alpha_colsum = arma::sum(alpha);
        alpha_sum += alpha_colsum;
        arma::fmat tmpx =  V_*alpha;
        X_sum += tmpx;
        arma::fmat alpha_2(alpha.n_rows,alpha.n_cols);
        for(int r=0;r<alpha_2.n_rows;++r)
        {
            arma::fmat tmpm = X_.each_col() - V_.col(r);
            alpha_2.row(r) = arma::sum(arma::square(tmpm));
        }
        arma::frowvec tmpvar = arma::sum(alpha_2%alpha);
        var_sum += tmpvar;
        ++idx;
    }
    X_ = X_sum.each_row() / alpha_sum;
    var = ( var_sum / (3.0*alpha_sum) );
    var += 1e-6;
    var = 1.0/var;
}


template<typename M>
void JRMPC<M>::stepMd()
{
    float mu = 0.0;
    arma::frowvec alpha_sumij(X_ptr->n_cols,arma::fill::zeros);
    std::vector<MatPtr>::iterator AIter;
    for(AIter=alpha_ptrs.begin();AIter!=alpha_ptrs.end();++AIter)
    {
        arma::fmat& alpha = **AIter;
        alpha_sumij += arma::sum(alpha);
        mu+=arma::sum(alpha_sumij);
    }
    mu*=(info_ptr->gamma+1.0);
    P_ = alpha_sumij.t();
    P_ /= mu;
}

template<typename M>
void JRMPC<M>::computeOnce(void)
{
    arma::fmat& X_ = *X_ptr;
    X_sum.fill(0.0);
    var_sum.fill(0.0);
    alpha_sum.fill(0.0);
    alpha_sumij.fill(0.0);
    T_updated_ = 0;
    X_updated_ = 0;
    float mu = 0.0;
    int idx = 0;
    std::vector<MatPtr>::iterator VIter;
    arma::fvec alpha_rowsum;
    arma::fmat U,V;
    arma::fvec s;
    arma::fmat C(3,3,arma::fill::eye);
    for(VIter=V_ptrs.begin();VIter!=V_ptrs.end();++VIter)
    {
        arma::fmat& V_ = **VIter;
        while(alpha_ptrs.size()<=idx)
        {
            alpha_ptrs.emplace_back(new arma::fmat(V_.n_cols,X_.n_cols));
        }
        //allocate alpha
        //compute alpha(step-E)
        arma::fmat& alpha = *alpha_ptrs[idx];
        #pragma omp parallel for
        for(int r=0;r<alpha.n_rows;++r)
        {
            alpha.row(r) = arma::sum(arma::square(X_.each_col() - V_.col(r)));
        }
        alpha.each_row()%=(-0.5*var);
        alpha = arma::trunc_exp(alpha);
        alpha.each_row()%=arma::pow(var,1.5)%P_;
        alpha_rowsum = arma::sum(alpha,1)+beta;
        alpha.each_col() /= alpha_rowsum;
        //update R t
        arma::fmat& R = *(res_ptr->Rs[idx]);
        arma::fmat dR;
        arma::fvec& t = *(res_ptr->ts[idx]);
        arma::fvec dt;
        arma::frowvec lambda = arma::sum( alpha );
        arma::fmat W = V_*alpha;
        W.each_row() %= var;
        arma::fvec mW = arma::sum(W,1);
        arma::frowvec b = lambda%var;
        arma::fvec mX = X_*b.t();
        float sumofWeights = arma::accu(b);
        arma::fmat A;
        A = X_*W.t() - mX*mW.t() / sumofWeights;
        dR = arma::fmat(3,3,arma::fill::eye);
        dt = arma::fvec(3,arma::fill::zeros);
        if(arma::svd(U,s,V,A,"std"))
        {
            C(2,2) = arma::det( U * V.t() )>=0 ? 1.0 : -1.0;
            dR = U*C*(V.t());
            dt = ( mX-dR*mW ) / sumofWeights;
        }
        R = dR*R;
        t = dR*t + dt;
        //evaluate if the T is really updated
        arma::fmat I(3,3,arma::fill::eye);
        if(
            ( info_ptr->eps_ < arma::accu(arma::square(dR - I)) ) ||
            ( info_ptr->eps_ < arma::accu(arma::square(dt)) )
           )
        {
            ++T_updated_;
        }

        //update V
        V_ = dR*V_;
        V_.each_col() += dt;

        //update X var
        alpha_sum += lambda;
        X_sum += V_*alpha;
        arma::fmat alpha_2(alpha.n_rows,alpha.n_cols);
        #pragma omp parallel for
        for(int r=0;r<alpha_2.n_rows;++r)
        {
            alpha_2.row(r) = arma::sum(arma::square(X_.each_col() - V_.col(r)));
        }
        arma::frowvec tmpvar = arma::sum(alpha_2%alpha);
        var_sum += tmpvar;
        ++idx;
    }
    //count how much X are updated
    arma::fmat newX = X_sum.each_row() / alpha_sum;
    arma::frowvec delta =  arma::sum( arma::square( X_ - newX ) );
    arma::uvec updatedX = arma::find( delta > info_ptr->eps_ );
    X_updated_ = updatedX.size();
    X_ = newX;
    var = ( (3.0*alpha_sum ) / ( var_sum + 3.0*alpha_sum*1e-6 ) );//restore reciprocal fractions of variation
    mu = arma::accu(alpha_sum);
    mu*=(info_ptr->gamma+1.0);
    P_ = alpha_sum;
    if( mu != 0)P_ /= mu;
}

template<typename M>
void JRMPC<M>::computeOnceJRCS(void)
{
    if(count==0)
    {
        int ridx =0;
        std::vector<std::shared_ptr<arma::fmat>>::iterator VIter;
        for( VIter = V_ptrs.begin() ; VIter != V_ptrs.end() ; ++VIter )
        {
            arma::fmat&v = **VIter;
            arma::fmat& R = *(res_ptr->Rs[ridx]);
            arma::fvec& t = *(res_ptr->ts[ridx]);
            v = R.i()*v;
            v.each_col() -= t;
            R.fill(arma::fill::eye);
            t*=-1.0;
            ++ridx;
        }
    }
    arma::fmat& X_ = *X_ptr;
    X_sum.fill(0.0);
    var_sum.fill(0.0);
    alpha_sum.fill(0.0);
    alpha_sumij.fill(0.0);
    T_updated_ = 0;
    X_updated_ = 0;
    float mu = 0.0;
    int idx = 0;
    std::vector<MatPtr>::iterator VIter;
    arma::fvec alpha_rowsum;
    arma::fmat U,V;
    arma::fvec s;
    arma::fmat C(3,3,arma::fill::eye);
    for(VIter=V_ptrs.begin();VIter!=V_ptrs.end();++VIter)
    {
        arma::fmat& V_ = **VIter;
        while(alpha_ptrs.size()<=idx)
        {
            alpha_ptrs.emplace_back(new arma::fmat(V_.n_cols,X_.n_cols));
        }
        //allocate alpha
        //compute alpha(step-E)
        arma::fmat& alpha = *alpha_ptrs[idx];
        arma::fmat& R = *(res_ptr->Rs[idx]);
        arma::fvec& t = *(res_ptr->ts[idx]);
        arma::fmat tX = R*X_;
        tX.each_col() += t;
        #pragma omp parallel for
        for(int r=0;r<alpha.n_rows;++r)
        {
            alpha.row(r) = arma::sum(arma::square(tX.each_col() - V_.col(r)));
        }
        alpha.each_row()%=(-0.5*var);
        alpha = arma::trunc_exp(alpha);
        alpha.each_row()%=arma::pow(var,1.5)%P_;
        alpha_rowsum = arma::sum(alpha,1)+beta;
        alpha.each_col() /= alpha_rowsum;
        //update R t
        arma::fmat dR;
        arma::fvec dt;
        arma::frowvec lambda = arma::sum( alpha );
        arma::fmat W = V_*alpha;
        W.each_row() %= var;
        arma::fvec mW = arma::sum(W,1);
        arma::frowvec b = lambda%var;
        arma::fvec mX = tX*b.t();
        float sumofWeights = arma::accu(b);
        arma::fmat A = W*tX.t() - mW*mX.t() / sumofWeights;
        dR = arma::fmat(3,3,arma::fill::eye);
        dt = arma::fvec(3,arma::fill::zeros);
        if(arma::svd(U,s,V,A,"std"))
        {
            C(2,2) = arma::det( U * V.t() )>=0 ? 1.0 : -1.0;
            dR = U*C*(V.t());
            dt = ( mW - dR*mX ) / sumofWeights;
        }
        R = dR*R;
        t = dR*t + dt;
        //evaluate if the T is really updated
        arma::fmat I(3,3,arma::fill::eye);
        if(
            ( info_ptr->eps_ < arma::accu(arma::square(dR - I)) ) ||
            ( info_ptr->eps_ < arma::accu(arma::square(dt)) )
           )
        {
            ++T_updated_;
        }

        //update V
        arma::fmat tV = V_.each_col() - t;
        tV = R.i()*tV;

        //update X var
        alpha_sum += lambda;
        X_sum += tV*alpha;
        arma::fmat alpha_2(alpha.n_rows,alpha.n_cols);
        #pragma omp parallel for
        for(int r=0;r<alpha_2.n_rows;++r)
        {
            alpha_2.row(r) = arma::sum(arma::square(X_.each_col() - tV.col(r)));
        }
        arma::frowvec tmpvar = arma::sum(alpha_2%alpha);
        var_sum += tmpvar;
        ++idx;
    }
    //count how much X are updated
    arma::fmat newX = X_sum.each_row() / alpha_sum;
    arma::frowvec delta =  arma::sum( arma::square( X_ - newX ) );
    arma::uvec updatedX = arma::find( delta > info_ptr->eps_ );
    X_updated_ = updatedX.size();
    X_ = newX;
    var = ( (3.0*alpha_sum ) / ( var_sum + 3.0*alpha_sum*1e-6 ) );//restore reciprocal fractions of variation
    mu = arma::accu(alpha_sum);
    mu*=(info_ptr->gamma+1.0);
    P_ = alpha_sum;
    if( mu != 0)P_ /= mu;
}

template<typename M>
void JRMPC<M>::setVarColor(uint32_t* color,int k)
{
    var_color = std::make_shared<arma::Col<uint32_t>>(color,k,false,true);
    varToColor();
}

template<typename M>
void JRMPC<M>::varToColor()
{
    if(var_color->size()==var.size())
    {
        float max_var = arma::max(var);
        float min_var = arma::min(var);
        float h;
        ColorArray::RGB32 tmp;
        int idx;
        for(idx=0;idx<var_color->size();++idx)
        {
            h = std::sqrt( ( var(idx) - min_var ) / ( max_var - min_var ) );
            ColorArray::hsv2rgb(h*220.0+5.0,0.5,1.0,tmp);
            (*var_color)(idx) = tmp.color;
        }
    }else{
        std::cerr<<"failed to trans color"<<std::endl;
    }
}

template<typename M>
void JRMPC<M>::restart(void)
{
    count = 0;
    arma::fmat rand_target = *X_ptr;
    initX(rand_target);
    *X_ptr = 0.9*(*X_ptr) + 0.1*rand_target;
}

template<typename M>
bool JRMPC<M>::isEnd()
{
    bool end = false;
    if( count >= info_ptr->max_iter )
    {
        if( restart_count >= info_ptr->max_restart )
        {
            std::cerr<<"N_Iter=="<<count<<std::endl;
            info_ptr->mode = MaxIter;
            end = true;
        }else{
            ++ restart_count;
            restart();
        }
    }
    if( 0==T_updated_ && count > 20)
    {
        std::cerr<<"T_updated_=="<<T_updated_<<std::endl;
        info_ptr->mode = Converge;
        end = true;
    }
    if( 0==X_updated_ && count > 20){
        std::cerr<<"X_updated_=="<<X_updated_<<std::endl;
        info_ptr->mode = Converge;
        end = true;
    }
    if( 0==X_updated_ && 0==T_updated_ && count > 20 )
    {
        std::cerr<<"X_updated_=="<<X_updated_<<"&&"<<"T_updated_=="<<T_updated_<<std::endl;
        info_ptr->mode = Converge;
        end = true;
    }
    if( end_ )
    {
        info_ptr->mode = Force;
        end = true;
    }
    if(end)
    {
        for(int idx=0;idx<res_ptr->Rs.size();++idx)
        {
            std::cout<< std::setprecision(8) << std::scientific;
            std::cout<<"R["<<idx<<"]=\n";
            arma::fmat& R = (*res_ptr->Rs[idx]);
            std::cout<<R(0,0)<<" "<<R(0,1)<<" "<<R(0,2)<<"\n";
            std::cout<<R(1,0)<<" "<<R(1,1)<<" "<<R(1,2)<<"\n";
            std::cout<<R(2,0)<<" "<<R(2,1)<<" "<<R(2,2)<<"\n";
            std::cout<<"\n";
        }
    }
    return end;
}

}
#endif

