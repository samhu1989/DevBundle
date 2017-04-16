#include "sjrcsbase.h"
#include "segmentationcore.h"
#include <cassert>
#include "iocore.h"
namespace JRCS{
SJRCSBase::SJRCSBase():JRCSBase()
{
    max_init_iter_ = 0;
    max_iter_ = 3;
    use_res_ = false;
}

bool SJRCSBase::configure(Config::Ptr config)
{
    config_ = config;
    if(config_->has("SJRCS_use_res"))
    {
        if(config_->getInt("SJRCS_use_res"))
        {
            use_res_ = true;
        }else{
            use_res_ = false;
        }
    }
    if(use_res_)
    {
        if(config_->has("SJRCS_res_step"))
        {
            res_step_ = config_->getInt("SJRCS_res_step");
        }else res_step_ = 2;
        if(config_->has("SJRCS_rebuild_k"))
        {
            rebuild_k_ = config_->getInt("SJRCS_rebuild_k");
        }else rebuild_k_ = 30;
        if(config_->has("SJRCS_base_k"))
        {
            base_k_ = config_->getInt("SJRCS_base_k");
        }else base_k_ = 90;
        if(config_->has("SJRCS_res_freq"))
        {
            res_act_freq_ = config_->getInt("SJRCS_res_freq");
        }else res_act_freq_ = 50;
    }
    return JRCSBase::configure(config);
}

bool SJRCSBase::isEnd(void)
{
    if(iter_count_>=max_iter_)return true;
    else return false;
}

void SJRCSBase::prepare_compute(void)
{
    iter_count_ = 1;
    reset_alpha();
    xtv_ptrlst_.resize(vvs_ptrlst_.size());
    #pragma omp for
    for( int i = 0 ; i < vvs_ptrlst_.size() ; ++i )
    {
        xtv_ptrlst_[i].reset(new arma::fmat(xv_ptr_->n_rows,xv_ptr_->n_cols,arma::fill::zeros));
    }
    reset_prob();
    if(use_res_)prepare_for_residue_correlation();
    QCoreApplication::processEvents();
    if(verbose_)std::cerr<<"done prepare"<<std::endl;
}

void SJRCSBase::step_a(int i)
{
    arma::fmat& xv_ = *xv_ptr_;
    arma::fmat& vv_ = *vvs_ptrlst_[i];
    //transform latent model
    arma::fmat& xtv_ = *xtv_ptrlst_[i];
    Ts& rt = rt_lst_[i];
    for(int o = 0 ; o < obj_num_ ; ++o )
    {
        arma::fmat R(rt[o].R,3,3,false,true);
        arma::fvec t(rt[o].t,3,false,true);
        arma::uword offset = 3*obj_range_[2*o];
        arma::uword size = obj_range_[2*o+1] - obj_range_[2*o] + 1;
        arma::fmat objv(((float*)xtv_.memptr())+offset,3,size,false,true);
        objv = R*xv_.cols(obj_range_[2*o],obj_range_[2*o+1]);
        objv.each_col() += t;
    }
    //calculate alpha
    arma::mat&  alpha = *alpha_ptrlst_[i];
    for(int r = 0 ; r < alpha.n_rows ; ++r )
    {
        arma::mat tmpv = arma::conv_to<arma::mat>::from(xtv_.each_col() - vv_.col(r));
        alpha.row(r)  = arma::sum(arma::square(tmpv));
    }
    alpha.each_row() %= (-0.5*x_invvar_);
    alpha = arma::trunc_exp(alpha);
    alpha.each_row() %= arma::pow(x_invvar_,1.5);
    alpha.each_row() %= x_p_;
    //normalise alpha
    arma::vec alpha_rowsum = ( 1.0 + beta_ ) * arma::sum(alpha,1);
    alpha.each_col() /= alpha_rowsum;

    //using residue correlation
    if( use_res_ && ( 0 == iter_count_% res_act_freq_) )
    {
        //truncating alpha with median to map the function
        arma::mat row_trunc_alpha;
        arma::vec trunc_alpha_rowsum;
        arma::mat col_trunc_alpha;
        arma::vec trunc_alpha_colsum;

        if(!funcs_.empty()) //only needed if we do calculate residue functions
        {
            truncate_alpha(
                        alpha,
                        row_trunc_alpha,
                        trunc_alpha_rowsum,
                        col_trunc_alpha,
                        trunc_alpha_colsum
                        );
        }

        //calculate residue functions
        for(int c = 0 ; c < funcs_.n_cols ; ++c )
        {
            arma::vec tmp;
            to_frame_func(row_trunc_alpha,trunc_alpha_rowsum,funcs_.col(c),tmp);
            to_vox_func(i,tmp,tmp);
            proj_and_rebuild(i,tmp,tmp);
            to_pix_func(i,tmp,tmp);
            to_model_func(col_trunc_alpha,trunc_alpha_colsum,tmp,tmp);
            res_(c,i) = res_energy(funcs_.col(c),tmp);
        }
    }
}

void SJRCSBase::step_b(void)
{
    if( use_res_ && ( 0 == iter_count_% res_act_freq_) )
    {
        if(!res_.empty())median_res_ = arma::mean(res_,1);
        if( verbose_ > 1 ){
            QString path;
            path = path.sprintf("./debug/res_cor/res_%u.mat",iter_count_);
            MATIO::save_to_matlab(res_,path.toStdString(),"res");
        }
        arma::rowvec res_err = arma::sum( arma::square( res_.each_col() - median_res_ ) );
        res_err.max(cir_frame_);
        cir_value_.fill(0);
     }
    QCoreApplication::processEvents();
}

void SJRCSBase::step_c(int i)
{
    arma::mat&  alpha = *alpha_ptrlst_[i];

    if( use_res_ && ( 0 == iter_count_% res_act_freq_ ) && ( i == cir_frame_ ) )
    {
        arma::uword cir = std::numeric_limits<arma::uword>::max();
        //max of circular correlation
        max_cir_cor(res_.col(i),median_res_,cir);
        //circulate the alpha accordingly
        cir_mat(cir,alpha);
        //recording cirvalue for debug
        cir_value_(i) = cir;
    }

    //update RT
    //#1 calculate weighted point cloud
    arma::fmat& vv = *vvs_ptrlst_[i];
    arma::fmat& vn = *vns_ptrlst_[i];
    arma::Mat<uint8_t>& vc = *vcs_ptrlst_[i];
    arma::fmat& wv = *wvs_ptrlst_[i];
    arma::fmat& wn = *wns_ptrlst_[i];
    arma::Mat<uint8_t>& wc = *wcs_ptrlst_[i];
    arma::frowvec alpha_colsum = arma::conv_to<arma::frowvec>::from(arma::sum( alpha ));
    calc_weighted(vv,vn,vc,alpha,wv,wn,wc);
    //#2 calculating RT for each object
    Ts& rt = rt_lst_[i];
    arma::fmat& xtv_ = *xtv_ptrlst_[i];
    for(int o = 0 ; o < obj_num_ ; ++o )
    {
        arma::fmat A;
        arma::fmat U,V;
        arma::fvec s;
        arma::fmat dR;
        arma::fvec dt;
        arma::fmat R(rt[o].R,3,3,false,true);
        arma::fvec t(rt[o].t,3,false,true);
        arma::uword offset = 3*obj_range_[2*o];
        arma::uword size = obj_range_[2*o+1] - obj_range_[2*o] + 1;
        arma::fmat objv(((float*)xtv_.memptr())+offset,3,size,false,true);
        arma::fmat v = wv.cols(obj_range_[2*o],obj_range_[2*o+1]);
        arma::fmat cv = v.each_col() - arma::mean(v,1);
        if(!v.is_finite())
        {
            std::cerr<<iter_count_<<":"<<o<<":!v.is_finite()"<<std::endl;
        }
        arma::frowvec square_lambda = arma::conv_to<arma::frowvec>::from(x_invvar_.cols(obj_range_[2*o],obj_range_[2*o+1]));
        square_lambda %= alpha_colsum.cols(obj_range_[2*o],obj_range_[2*o+1]);
        square_lambda += std::numeric_limits<float>::epsilon(); // for numeric stability
        if(!square_lambda.is_finite())
        {
            std::cerr<<"!square_lambda.is_finite()"<<std::endl;
        }
        arma::frowvec p(square_lambda.n_cols,arma::fill::ones);
        arma::frowvec square_norm_lambda = square_lambda / arma::accu(square_lambda);
        if(!square_norm_lambda.is_finite())
        {
            std::cerr<<"!square_norm_lambda.is_finite()"<<std::endl;
        }
        p -= square_norm_lambda;
        arma::fmat tmp = objv.each_col() - arma::mean(objv,1) ;
        tmp.each_row() %= p % square_lambda;
        A = cv*tmp.t();
        switch(rttype_)
        {
        case Gamma:
        {
            arma::fmat B = A.submat(0,0,1,1);
            dR = arma::fmat(3,3,arma::fill::eye);
            if(arma::svd(U,s,V,B,"std"))
            {
                arma::fmat C(2,2,arma::fill::eye);
                C(1,1) = arma::det( U * V.t() )>=0 ? 1.0 : -1.0;
                arma::fmat dR2D = U*C*(V.t());
                dR.submat(0,0,1,1) = dR2D;
                if(!dR.is_finite())
                {
                    std::cerr<<iter_count_<<":!dR.is_finite()"<<std::endl;
                    dR.fill(arma::fill::eye);
                }
                arma::fmat tmp = v - dR*objv;
                tmp.each_row() %= square_norm_lambda;
                dt = arma::sum(tmp,1);
                if(!dt.is_finite())
                {
                    std::cerr<<iter_count_<<":!dt.is_finite()"<<std::endl;
                    dt.fill(arma::fill::zeros);
                }
            }
        }
            break;
        default:
        {
            if(arma::svd(U,s,V,A,"std"))
            {
                arma::fmat C(3,3,arma::fill::eye);
                C(2,2) = arma::det( U * V.t() )>=0 ? 1.0 : -1.0;
                dR = U*C*(V.t());
                if(!dR.is_finite())
                {
                    std::cerr<<iter_count_<<":!dR.is_finite()"<<std::endl;
                    dR.fill(arma::fill::eye);
                }
                arma::fmat tmp = v - dR*objv;
                tmp.each_row() %= square_norm_lambda;
                dt = arma::sum(tmp,1);
                if(!dt.is_finite())
                {
                    std::cerr<<iter_count_<<":!dt.is_finite()"<<std::endl;
                    dt.fill(arma::fill::zeros);
                }
            }
        }
        }

        objv = dR*objv;
        objv.each_col() += dt;

        //updating R T
        R = dR*R;
        t = dR*t + dt;
    }

    arma::fmat alpha_2(alpha.n_rows,alpha.n_cols);
    for(int r=0;r<alpha_2.n_rows;++r)
    {
        alpha_2.row(r) = arma::sum(arma::square(xtv_.each_col() - vv.col(r)));
    }
    vvar_.row(i) = arma::sum(alpha_2%alpha);
//    std::cerr<<"step_c("<<i<<")"<<std::endl;
}

void SJRCSBase::step_d(void)
{
    if( use_res_ && ( 0 == iter_count_% res_act_freq_) )
    {
        if( verbose_ > 1 ){
            QString path;
            path = path.sprintf("./debug/res_cor/cir_%u.mat",iter_count_);
            MATIO::save_to_matlab(cir_value_,path.toStdString(),"cir");
        }
        reset_x();
        reset_var_p();
        return;
    }
    //Updating X
    xv_ptr_->fill(0.0);
    xn_ptr_->fill(0.0);
    arma::fmat xc(xc_ptr_->n_rows,xc_ptr_->n_cols,arma::fill::zeros);
    for(int i=0;i<wvs_ptrlst_.size();++i)
    {
        float rate0 = float(i) / float(i+1.0);
        float rate1 = 1.0 - rate0;
        if( i >= rt_lst_.size() )
        {
            std::cerr<<"i >= rt_lst_.size()"<<std::endl;
        }
        Ts& rt = rt_lst_[i];
        arma::fmat& wv = *wvs_ptrlst_[i];
        arma::fmat& wn = *wns_ptrlst_[i];
        arma::fmat wc = arma::conv_to<arma::fmat>::from(*wcs_ptrlst_[i]);
        #pragma omp parallel for
        for(int o = 0 ; o < obj_num_ ; ++o )
        {
            arma::fmat R(rt[o].R,3,3,false,true);
            arma::fvec t(rt[o].t,3,false,true);
            arma::fmat invR = R.i();
            if(!invR.is_finite())
            {
                std::cerr<<iter_count_<<":"<<o<<":!invR.is_finite()"<<std::endl;
            }
            if(!t.is_finite())
            {
                std::cerr<<iter_count_<<":"<<o<<":!t.is_finite()"<<std::endl;
            }
            arma::uword s = obj_range_[2*o];
            arma::uword e = obj_range_[2*o+1];
            arma::fmat tv = wv.cols(s,e);
            if(!tv.is_finite())
            {
                std::cerr<<iter_count_<<":"<<o<<":!tv.is_finite()"<<std::endl;
            }
            tv.each_col() -= t;
            xv_ptr_->cols(s,e) =  rate0*xv_ptr_->cols(s,e)+rate1*(invR*tv);
            xn_ptr_->cols(s,e) = rate0*xv_ptr_->cols(s,e)+rate1*invR*wn.cols(s,e);
            xc.cols(s,e) = rate0*xc.cols(s,e)+rate1*wc.cols(s,e);
        }
    }
    if(!xv_ptr_->is_finite())
    {
        std::cerr<<iter_count_<<":!xv_ptr_->is_finite()"<<std::endl;
    }
    if(!xn_ptr_->is_finite())
    {
        std::cerr<<iter_count_<<":!xn_ptr_->is_finite()"<<std::endl;
    }
    if(!xc.is_finite())
    {
        std::cerr<<"!xc.is_finite()"<<std::endl;
    }
    //fix the x center position
//    std::cerr<<"xv:"<<xv_ptr_->n_cols<<std::endl;
    #pragma omp parallel for
    for(int o = 0 ; o < obj_num_ ; ++o )
    {
        arma::uword s = obj_range_[2*o];
        arma::uword e = obj_range_[2*o+1];
        arma::fmat newxv = xv_ptr_->cols(s,e);
        arma::fvec t =  obj_pos_.col(o) - arma::mean(newxv,1);
        xv_ptr_->cols(s,e) = newxv.each_col() + t;
    }
    *xn_ptr_ = arma::normalise(*xn_ptr_);
    *xc_ptr_ = arma::conv_to<arma::Mat<uint8_t>>::from( xc );

    //Updating variance of centroids
    arma::rowvec alpha_sum;
    for(int i=0;i<alpha_ptrlst_.size();++i)
    {
        arma::mat& alpha = *alpha_ptrlst_[i];
        if(alpha_sum.empty())alpha_sum = arma::sum(alpha);
        else alpha_sum += arma::sum(alpha);
    }
    x_invvar_ = ( ( 3.0*alpha_sum ) / ( arma::sum(vvar_) + beta_ ) );//restore reciprocal fractions of variation
    //Updating weight of centroids
    float mu = arma::accu(alpha_sum);
    mu *= ( 1.0 + beta_ );
    x_p_ = alpha_sum;
    if( mu != 0)x_p_ /= mu;
}

void SJRCSBase::finish_steps(void)
{
    if(verbose_>1)std::cerr<<"finishing steps"<<std::endl;
    update_color_label();
    if(verbose_>1)std::cerr<<"color label updated"<<std::endl;
    calc_obj();
    std::cout<<"obj("<<iter_count_<<"):"<<-0.5*obj_vec_.back()<<std::endl;
    ++ iter_count_;
    QCoreApplication::processEvents();
}

//prepare circulated indicating functions
void SJRCSBase::prepare_for_residue_correlation(void)
{
    std::cerr<<"preparing residue correlation"<<std::endl;
    arma::vec func_(xv_ptr_->n_cols,arma::fill::zeros);
    #pragma omp parallel for
    for( int o = 0 ; o < obj_num_ ; ++ o)
    {
        arma::uword obj_size = obj_range_[2*o+1] - obj_range_[2*o] + 1;
        func_.rows(obj_range_[2*o],obj_range_[2*o+1]) = ( 1.0 + arma::cos(arma::linspace<arma::vec>(-M_PI,M_PI,obj_size)) ) / 2.0 / double (o + 1);
    }
    arma::uvec cir_indice_ = arma::linspace<arma::uvec>(0,xv_ptr_->n_cols-1,xv_ptr_->n_cols);
    arma::uword N = func_.size();
    arma::uword N_step = std::floor( xv_ptr_->n_cols / res_step_ ) ;
    cir_pos_ = arma::linspace<arma::uvec>( 0 , func_.size() - 1 , N_step );
    funcs_ = arma::mat(func_.size(),cir_pos_.size(),arma::fill::zeros);
    cir_indices_ = arma::umat(func_.size(),cir_pos_.size(),arma::fill::zeros);
    std::cerr<<"a"<<std::endl;
    #pragma omp parallel for
    for( int i = 0 ; i < N_step ; ++i )
    {
        funcs_.submat(cir_pos_(i),i,N-1,i) = func_.rows(0,N-1-cir_pos_(i));
        cir_indices_.submat(cir_pos_(i),i,N-1,i) = cir_indice_.rows(0,N-1-cir_pos_(i));
        if( cir_pos_(i) > 0 ){
            funcs_.submat(0,i,cir_pos_(i)-1,i) = func_.rows( N - cir_pos_(i) , N - 1 );
            cir_indices_.submat(0,i,cir_pos_(i)-1,i) = cir_indice_.rows( N - cir_pos_(i) , N - 1 );
        }
    }
    std::cerr<<"b"<<std::endl;
    cir_value_ = arma::uvec(vvs_ptrlst_.size(),arma::fill::zeros);
    res_ = arma::mat(cir_pos_.size(),vvs_ptrlst_.size(),arma::fill::zeros);
}

void SJRCSBase::truncate_alpha(
        const arma::mat& alpha,
        arma::mat& row_trunc_alpha,
        arma::vec& trunc_alpha_rowsum,
        arma::mat& col_trunc_alpha,
        arma::vec& trunc_alpha_colsum
        )
{

    //row_trunc_alpha for forward mapping
    arma::vec alpha_row_median = arma::max( alpha , 1 );
    row_trunc_alpha = alpha;
    for(int r=0 ; r < row_trunc_alpha.n_rows ; ++r )
    {
        arma::rowvec row = row_trunc_alpha.row(r);
        row( row < alpha_row_median(r) ).fill(0.0);
        row_trunc_alpha.row(r) = row;
    }
    row_trunc_alpha += std::numeric_limits<double>::epsilon(); //add eps for numeric stability
    trunc_alpha_rowsum = arma::sum( row_trunc_alpha , 1 );

    //col_trunc_alpha for backward mapping
    arma::rowvec alpha_col_median = arma::max( alpha );
    col_trunc_alpha = alpha;
    for(int c=0 ; c < col_trunc_alpha.n_cols ; ++c )
    {
        arma::vec col = col_trunc_alpha.col(c);
        col( col < alpha_col_median(c) ).fill(0.0);
        col_trunc_alpha.col(c) = col;
    }
    col_trunc_alpha += std::numeric_limits<double>::epsilon(); //add eps for numeric stability
    trunc_alpha_colsum = arma::sum( col_trunc_alpha ).t();
}

void SJRCSBase::to_frame_func(
        const arma::mat& row_trunc_alpha,
        const arma::vec& trunc_alpha_rowsum,
        const arma::vec& model_func,
        arma::vec& frame_func)
{
    if(model_func.n_rows!=row_trunc_alpha.n_cols)return;
    frame_func = row_trunc_alpha*model_func;
    for(int r=0;r<row_trunc_alpha.n_rows;++r)
    {
        if( 0 != trunc_alpha_rowsum(r) )
        {
            frame_func(r) /= trunc_alpha_rowsum(r);
        }
    }
}

void SJRCSBase::to_vox_func(int index,const arma::vec& frame_func,arma::vec& vox_func)
{
    MeshBundle<DefaultMesh>::Ptr m_ptr = inputs_[index];
    arma::vec tmp = frame_func;
    m_ptr->graph_.getVoxFunc(tmp,vox_func);
}

void SJRCSBase::proj_and_rebuild(int index,const arma::vec& vox_func,arma::vec& re_vox_func)
{
    arma::vec coeff = arma::vectorise(vox_func.t()*( *bases_[index] ));
    arma::uvec sorted_i = arma::sort_index(arma::abs(coeff),"descend");
    arma::uvec selected_i = sorted_i.head(rebuild_k_);
    arma::vec sub_coeff = coeff.rows(selected_i);
    arma::mat sub_bases = (*bases_[index]).cols(selected_i);
    re_vox_func = sub_bases*sub_coeff;
}

void SJRCSBase::to_pix_func(int index,const arma::vec& vox_func,arma::vec& pix_func)
{
    MeshBundle<DefaultMesh>::Ptr m_ptr = inputs_[index];
    arma::vec tmp = vox_func;
    m_ptr->graph_.getPixFunc(tmp,pix_func);
}

void SJRCSBase::to_model_func(
        const arma::mat& col_trunc_alpha,
        const arma::vec& trunc_alpha_colsum,
        const arma::vec& frame_func,
        arma::vec& model_func)
{
    if(frame_func.n_rows!=col_trunc_alpha.n_rows)return;
    model_func = (col_trunc_alpha.t())*frame_func;
    for(int c=0;c<col_trunc_alpha.n_cols;++c)
    {
        if( 0 != trunc_alpha_colsum(c) )
        {
            model_func(c) /= trunc_alpha_colsum(c);
        }
    }
}

double SJRCSBase::res_energy(const arma::vec& func0,const arma::vec& func1)
{
    if(func0.size()!=func1.size()){
        std::cerr<<"func0.size()!=func1.size()"<<std::endl;
        return -1.0;
    }
    return arma::accu(arma::square(func0 - func1));
}

void SJRCSBase::max_cir_cor(const arma::vec& func0,const arma::vec& func1,arma::uword& max_cir)
{
    if(func0.size()!=func1.size()){
        std::cerr<<"func0.size()!=func1.size()"<<std::endl;
        return;
    }
    arma::vec cor(func0.size(),arma::fill::zeros);
    arma::uword N = func0.size();
    double max = std::numeric_limits<double>::lowest();
    for(arma::uword i=0 ; i<N ; ++i)//position of func0
        for(arma::uword j=0 ; j<N ; ++j)//circulation
        {
            arma::uword k = (i+j) % N ;
            cor(j) += func0(i)*func1(k);
            if( cor(j) > max ){
                max = cor(j);
                max_cir = j;
            }
        }
}

void SJRCSBase::cir_mat(const arma::uword& cir,arma::mat& alpha)
{
    //reverse the circulation
    alpha.cols(cir_indices_.col(cir)) = alpha;
}

void SJRCSBase::calc_weighted(const arma::fmat&vv,
        const arma::fmat&vn,
        arma::Mat<uint8_t>&vc,
        const arma::mat &alpha,
        arma::fmat&wv,
        arma::fmat&wn,
        arma::Mat<uint8_t>&wc
        )
{
    arma::fmat fwc = arma::conv_to<arma::fmat>::from(wc);
    arma::rowvec alpha_median = arma::median( alpha );

    arma::mat trunc_alpha = alpha;
    for(int c=0;c<alpha.n_cols;++c)
    {
        arma::vec col = trunc_alpha.col(c);
        col( col < alpha_median(c) ).fill(0.0);
        trunc_alpha.col(c) = col;
    }

    trunc_alpha += std::numeric_limits<double>::epsilon(); //add eps for numeric stability

    arma::rowvec trunc_alpha_colsum = arma::sum( trunc_alpha );

    wv = vv*arma::conv_to<arma::fmat>::from(trunc_alpha);

    for(int c=0;c<wn.n_cols;++c)
    {
//        std::cerr<<"wn("<<c<<")"<<std::endl;
        arma::mat wvn = arma::conv_to<arma::mat>::from(vn);
        wvn.each_row() %= trunc_alpha.col(c).t();
        arma::vec eigval;
        arma::mat eigvec;
        arma::eig_sym(eigval,eigvec,wvn*wvn.t(),"std");
        wn.col(c) = arma::conv_to<arma::fvec>::from(eigvec.col(2));
//        std::cerr<<wn.col(c)<<std::endl;
    }

    fwc = arma::conv_to<arma::fmat>::from(vc)*arma::conv_to<arma::fmat>::from(trunc_alpha);

    if(!wv.is_finite())
    {
        std::cerr<<iter_count_<<":a:!wv.is_finite()"<<std::endl;
    }

    for(int c=0;c<alpha.n_cols;++c)
    {
        if( 0 != trunc_alpha_colsum(c) )
        {
            wv.col(c) /= trunc_alpha_colsum(c);
//            wn.col(c) /= trunc_alpha_colsum(c);
            fwc.col(c) /= trunc_alpha_colsum(c);
        }
    }

    if(!wv.is_finite())
    {
        std::cerr<<iter_count_<<":b:!wv.is_finite()"<<std::endl;
    }

//    wn = arma::normalise( wn );
    assert(wv.is_finite());
    assert(wn.is_finite());
    assert(fwc.is_finite());
    wc = arma::conv_to<arma::Mat<uint8_t>>::from(fwc);
}

int SJRCSBase::evaluate_k()
{
    std::cerr<<"SJRCSBase::evaluate_k()"<<std::endl;
    if(vvs_ptrlst_.empty())
    {
        throw std::logic_error("need input v before evaluate_k");
    }
    MatPtrLst::iterator iter;
    arma::uvec k_lst(vvs_ptrlst_.size());
    int idx = 0;
    for(iter=vvs_ptrlst_.begin();iter!=vvs_ptrlst_.end();++iter)
    {
        k_lst(idx) = (*iter)->n_cols;
        ++idx;
    }
    arma::uword k = arma::median(k_lst);
    k = k/2 + 7 ; //based on median size but at least 7;
    std::cerr<<"SJRCSBase::evaluate_k():k="<<int(k)<<std::endl;
    if(k<0)
    {
        std::cerr<<"k_lst:"<<k_lst.t()<<std::endl;
    }
    return int(k);
}

void SJRCSBase::initx(
        const MatPtr& xv,
        const MatPtr& xn,
        const CMatPtr& xc
        )
{
    if(verbose_)std::cerr<<"SJRCSBase start init x"<<std::endl;
    xv_ptr_ = xv;
    xn_ptr_ = xn;
    xc_ptr_ = xc;

    max_obj_radius_ = 1.0;// for randomly reset rt
    init_alpha_ = false;  //for this method no init of alpha is used

    int k = xv_ptr_->n_cols;
    if(verbose_){
        std::cerr<<"k:"<<k<<"=="<<xn_ptr_->n_cols<<"=="<<xc_ptr_->n_cols<<std::endl;
    }
    int r_k = k;
    float* pxv = (float*)xv_ptr_->memptr();
    float* pxn = (float*)xn_ptr_->memptr();
    uint8_t* pxc = (uint8_t*)xc_ptr_->memptr();
    int N = obj_prob_.size();
    obj_pos_ = arma::fmat(3,N,arma::fill::zeros);
    arma::frowvec z = arma::linspace<arma::frowvec>(float(-N/2),float(N/2),N);
    obj_pos_.row(2) = z;
    if(verbose_>0){
        std::cerr<<"obj_pos_:"<<std::endl;
        std::cerr<<obj_pos_<<std::endl;
    }
    obj_range_.resize(2*N);
    arma::uword* pxr = (arma::uword*)obj_range_.data();
    arma::uword* pxr_s = pxr;
    *pxr = 0;
    for(int obj_idx = 0 ; obj_idx < obj_prob_.size() ; ++ obj_idx )
    {
        int obj_size = int(float(k)*float(obj_prob_(obj_idx)));
        obj_size = std::max(9,obj_size);
        obj_size = std::min(r_k,obj_size);
        if( obj_idx == (obj_prob_.size() - 1 ) ) obj_size = r_k;

        if( pxr > pxr_s )*pxr = (*(pxr - 1))+1;
        *(pxr+1) = *pxr + obj_size - 1;

        arma::fmat objv(pxv,3,obj_size,false,true);
        arma::fmat objn(pxn,3,obj_size,false,true);
        arma::Mat<uint8_t> objc(pxc,3,obj_size,false,true);

        pxv += 3*obj_size;
        pxn += 3*obj_size;
        pxc += 3*obj_size;
        pxr += 2;
        r_k -= obj_size;

        arma::fvec pos = obj_pos_.col(obj_idx);
        reset_obj_vn(0.5,pos,objv,objn);
        reset_obj_c(objc);
    }
    if(verbose_)std::cerr<<"Done initx_"<<std::endl;
}

void SJRCSBase::reset_x(void)
{
    float* pxv = (float*)xv_ptr_->memptr();
    float* pxn = (float*)xn_ptr_->memptr();
    uint8_t* pxc = (uint8_t*)xc_ptr_->memptr();
    for(int o = 0 ; o < obj_prob_.size() ; ++ o )
    {
        arma::uword obj_size = obj_range_[2*o+1] - obj_range_[2*o] + 1;
        arma::fmat objv(pxv,3,obj_size,false,true);
        arma::fmat objn(pxn,3,obj_size,false,true);
        arma::Mat<uint8_t> objc(pxc,3,obj_size,false,true);

        arma::fmat v = xv_ptr_->cols(obj_range_[2*o],obj_range_[2*o]+1);
        arma::fvec maxVXYZ = arma::max(v,1);
        arma::fvec minVXYZ = arma::min(v,1);

        float r = maxVXYZ(2) - minVXYZ(2);
        r = std::sqrt(r*r);

        pxv += 3*obj_size;
        pxn += 3*obj_size;
        pxc += 3*obj_size;

        arma::fvec pos = obj_pos_.col(o);
        reset_obj_vn(std::max(r,0.1f),pos,objv,objn);
        reset_obj_c(objc);
    }
}

void SJRCSBase::rand_sphere(
        arma::fmat& ov
        )
{
    int k = ov.n_cols;
    arma::frowvec az = arma::randu<arma::frowvec>(k);
    arma::frowvec el = arma::randu<arma::frowvec>(k);
    az*=2*M_PI;
    el*=2*M_PI;
    ov.row(0) = arma::cos(az)%arma::cos(el);
    ov.row(1) = arma::sin(el);
    ov.row(2) = arma::sin(az)%arma::cos(el);
    #pragma omp parallel for
    for(int i=0 ; i < ov.n_cols ; ++i )
    {
        std::random_shuffle(ov.begin_col(i),ov.end_col(i));
    }
}

void SJRCSBase::reset_prob()
{
    if(verbose_>0)std::cerr<<"initializing probability"<<std::endl;

    arma::fvec maxAllXYZ,minAllXYZ;
    MatPtrLst::iterator vviter;
    for( vviter = vvs_ptrlst_.begin() ; vviter != vvs_ptrlst_.end() ; ++vviter )
    {
        arma::fmat&v = **vviter;
        arma::fvec maxVXYZ = arma::max(v,1);
        arma::fvec minVXYZ = arma::min(v,1);
        if(vviter==vvs_ptrlst_.begin())
        {
            maxAllXYZ = maxVXYZ;
            minAllXYZ = minVXYZ;
        }else{
            maxAllXYZ = arma::max(maxVXYZ,maxAllXYZ);
            minAllXYZ = arma::min(minVXYZ,minAllXYZ);
        }
    }

    float r = maxAllXYZ(2)-minAllXYZ(2);
    float maxvar = r*r;
    x_invvar_ = arma::rowvec(xv_ptr_->n_cols);
    x_invvar_.fill(1.0/maxvar);

    vvar_ = arma::mat(vvs_ptrlst_.size(),xv_ptr_->n_cols);

    x_p_ = arma::rowvec(xv_ptr_->n_cols);
    x_p_.fill(1.0/float(xv_ptr_->n_cols));

    if(verbose_>0)std::cerr<<"done probability"<<std::endl;
}

void SJRCSBase::reset_var_p(void)
{
    for(int o=0;o<obj_prob_.size();++o)
    {
        arma::fmat v = xv_ptr_->cols(obj_range_[2*o],obj_range_[2*o]+1);
        arma::fvec maxVXYZ = arma::max(v,1);
        arma::fvec minVXYZ = arma::min(v,1);
        float maxvar = arma::accu(arma::square(maxVXYZ-minVXYZ));
        x_invvar_.cols(obj_range_[2*o],obj_range_[2*o]+1).fill(1.0/maxvar);
        arma::rowvec po = x_p_.cols(obj_range_[2*o],obj_range_[2*o]+1);
        x_p_.cols(obj_range_[2*o],obj_range_[2*o]+1).fill(arma::mean(po));
    }
}

void SJRCSBase::update_color_label()
{
    for(int idx=0;idx<vvs_ptrlst_.size();++idx)
    {
        arma::mat& alpha = *alpha_ptrlst_[idx];
        arma::mat obj_p(alpha.n_rows,obj_num_);
        arma::Col<uint32_t>& vl = *vls_ptrlst_[idx];
        #pragma omp parallel for
        for(int o = 0 ; o < obj_num_ ; ++o )
        {
            arma::mat sub_alpha = alpha.cols(obj_range_[2*o],obj_range_[2*o+1]);
            obj_p.col(o) = arma::sum(sub_alpha,1);
        }
        arma::uvec label(alpha.n_rows);
        #pragma omp parallel for
        for(int r = 0 ; r < obj_p.n_rows ; ++r )
        {
            arma::uword l;
            arma::rowvec point_prob = obj_p.row(r);
            point_prob.max(l);
            label(r) = l+1;
        }
        ColorArray::colorfromlabel((uint32_t*)vl.memptr(),vl.size(),label);
    }
}

void SJRCSBase::calc_obj(void)
{
    if(verbose_>1)std::cerr<<"SJRCSBase::calc_obj()"<<std::endl;
    obj_vec_.push_back(0.0);
    for(int idx=0;idx<vvs_ptrlst_.size();++idx)
    {
        arma::mat& alpha = *alpha_ptrlst_[idx];
        arma::fmat& vv_ = *vvs_ptrlst_[idx];
        Ts& rt = rt_lst_[idx];
        arma::fmat& xtv_ = *xtv_ptrlst_[idx];
        arma::fmat alpha_2(alpha.n_rows,alpha.n_cols);
        #pragma omp parallel for
        for(int r=0;r<alpha_2.n_rows;++r)
        {
            alpha_2.row(r) = arma::sum(arma::square(xtv_.each_col() - vv_.col(r)));
            alpha_2.row(r) %= arma::conv_to<arma::frowvec>::from(x_invvar_);
            alpha_2.row(r) -= 1.5*arma::conv_to<arma::frowvec>::from(arma::trunc_log(x_invvar_));
            alpha_2.row(r) -= 2.0*arma::conv_to<arma::frowvec>::from(arma::trunc_log(x_p_));
        }
        obj_vec_.back()+=arma::accu(alpha%alpha_2);
    }
    obj_vec_.back()*=-0.5;
}

bool SJRCSBase::input_extra(
        const MeshBundle<DefaultMesh>::PtrList& inputs
        )
{
    if(!use_res_)return true;
    std::cerr<<"input_extra:"<<std::endl;
    inputs_ = inputs;
    if(inputs_.empty())return false;
    if(inputs_[0]->graph_.empty())return false;
    bases_.resize(inputs_.size());
//    #pragma omp parallel for
    for( int i=0 ; i < inputs_.size() ; ++i )
    {
        MeshBundle<DefaultMesh>::Ptr m_ptr = inputs_[i];
        Segmentation::NormalizedCuts<DefaultMesh> ncuts_;
        assert(ncuts_.configure(config_));
        ncuts_.setK( base_k_ + 5 );
        ncuts_.setEpsilon( std::numeric_limits<double>::lowest() );
        ncuts_.setType(Segmentation::NormalizedCuts<DefaultMesh>::M);
        ncuts_.computeW_Graph(m_ptr);
        ncuts_.decompose();
        bases_[i].reset(new arma::mat());
        (*bases_[i]) = ncuts_.getY().cols(0,base_k_);
    }
    std::cerr<<"end input_extra:"<<std::endl;
}

}
