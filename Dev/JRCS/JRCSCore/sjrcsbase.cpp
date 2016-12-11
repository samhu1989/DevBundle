#include "sjrcsbase.h"
namespace JRCS{
SJRCSBase::SJRCSBase():max_iter_(30),beta_(std::numeric_limits<double>::epsilon())
{
    arma::arma_rng::set_seed_random();
}
bool SJRCSBase::configure(Config::Ptr config)
{
    config_ = config;
}
void SJRCSBase::input(
        const MatPtrLst& vv,
        const MatPtrLst& vn,
        const CMatPtrLst& vc,
        const LCMatPtrLst& vl
        )
{
    ;
}
void SJRCSBase::resetw(
        const MatPtrLst& wv,
        const MatPtrLst& wn,
        const CMatPtrLst& wc
        )
{
    ;
}
void SJRCSBase::initx(
        const MatPtr& xv,
        const MatPtr& xn,
        const CMatPtr& xc
        )
{
    ;
}
bool SJRCSBase::isEnd(void)
{
    if(iter_count_>max_iter_)return true;
}
void SJRCSBase::prepare_compute(void)
{
    iter_count_ = 0;
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
        arma::uword offset = 3*obj_range_[2*o]+1;
        arma::uword size = obj_range_[2*o+1] - obj_range_[2*o] + 1;
        arma::fmat objv(((float*)xtv_.memptr())+offset,3,size,false,true);
        objv = R*xv_.cols(obj_range_[2*o],obj_range_[2*o]+1);
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
    alpha = arma::exp(alpha);
    alpha.each_row() %= arma::pow(x_invvar_,1.5);
    alpha.each_row() %= x_p_;

    //normalise alpha
    arma::vec alpha_rowsum = ( 1.0 + beta_ ) * arma::sum(alpha,1);
    alpha.each_col() /= alpha_rowsum;
    alpha_rowsum = arma::sum(alpha,1);

    //calculate residue functions
    for(int c = 0 ; c < funcs_.n_cols ; ++c )
    {
        arma::vec tmp;
        to_frame_func(alpha,funcs_.col(c),tmp);
        to_vox_func(i,tmp,tmp);
        proj_and_rebuild(i,tmp,tmp);
        to_pix_func(i,tmp,tmp);
        to_model_func(alpha,tmp,tmp);
        res_(c,i) = res_energy(funcs_.col(c),tmp);
    }
}

void SJRCSBase::step_b(void)
{
    median_res_ = arma::median(res_,1);
}

void SJRCSBase::step_c(int i)
{
    //max of circular correlation
    arma::uword cir;
    max_cir_cor(median_res_,res_.col(i),cir);
    arma::mat&  alpha = *alpha_ptrlst_[i];
    //circulate the alpha accordingly
    cir_mat(cir,alpha);

    //update RT
    //#1 calculate weighted point cloud
    arma::fmat& vv = *vvs_ptrlst_[i];
    arma::fmat& vn = *vns_ptrlst_[i];
    arma::Mat<uint8_t>& vc = *vcs_ptrlst_[i];
    arma::fmat& wv = *wvs_ptrlst_[i];
    arma::fmat& wn = *wns_ptrlst_[i];
    arma::Mat<uint8_t>& wc = *wcs_ptrlst_[i];
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
        arma::uword offset = 3*obj_range_[2*o]+1;
        arma::uword size = obj_range_[2*o+1] - obj_range_[2*o] + 1;
        arma::fmat objv(((float*)xtv_.memptr())+offset,3,size,false,true);
        arma::fmat v = wv.cols(obj_range_[2*o],obj_range_[2*o+1]);
        arma::fmat cv = v.each_col() - arma::mean(v,1);
        A = cv*( objv.each_col() - arma::mean(objv,1) );
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
                dt = arma::mean( v - dR*objv,1);
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
                dt = arma::mean( v - dR*objv,1);
            }
        }
        }

        //updating R T
        R = dR*R;
        t = dR*t + dt;
    }

    arma::fmat alpha_2(alpha.n_rows,alpha.n_cols);
    for(int r=0;r<alpha_2.n_rows;++r)
    {
        alpha_2.row(r) = arma::sum(arma::square(xtv_.each_col() - vv.col(r)));
    }
    var_.row(i) = arma::sum(alpha_2%alpha);

}
void SJRCSBase::step_d(void)
{
    //Updating X
    xv_ptr_->fill(0.0);
    xn_ptr_->fill(0.0);
    arma::fmat xc(xc_ptr_->n_rows,xc_ptr_->n_cols,arma::fill::zeros);
    for(int i=0;i<wvs_ptrlst_.size();++i)
    {
        float rate0 = float(i) / float(i+1);
        float rate1 = 1.0 - rate0;
        Ts& rt = rt_lst_[i];
        arma::fmat& wv = *wvs_ptrlst_[i];
        arma::fmat& wn = *wns_ptrlst_[i];
        arma::fmat wc = arma::conv_to<arma::fmat>::from(*wcs_ptrlst_[i]);
        #pragma omp parallel for
        for(int o = 0 ; o < obj_num_ ; ++o )
        {
            arma::fmat R(rt[o].R,3,3,false,true);
            arma::fvec t(rt[o].t,3,false,true);
            arma::uword s = obj_range_[2*o];
            arma::uword e = obj_range_[2*o+1];
            arma::fmat tv = wv.cols(s,e);
            tv.each_col() -= t;
            xv_ptr_->cols(s,e) =  rate0*xv_ptr_->cols(s,e)+rate1*(R.i()*tv);
            xn_ptr_->cols(s,e) = rate0*xv_ptr_->cols(s,e)+rate1*R.i()*wn.cols(s,e);
            xc.cols(s,e) = rate0*xc.cols(s,e)+rate1*wc.cols(s,e);
        }
    }
    //fix the x center position
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
    x_invvar_ = ( ( 3.0*alpha_sum ) / ( arma::sum(var_) + beta_ ) );//restore reciprocal fractions of variation
    //Updating weight of centroids
    float mu = arma::accu(alpha_sum);
    mu *= ( 1.0 + beta_ );
    x_p_ = alpha_sum;
    if( mu != 0)x_p_ /= mu;
}
void SJRCSBase::finish_steps(void)
{
    ++ iter_count_;
}
void SJRCSBase::to_frame_func(const arma::mat& alpha,const arma::vec& model_func,arma::vec& frame_func)
{
    ;
}
void SJRCSBase::to_vox_func(int index,const arma::vec& frame_func,arma::vec& vox_func)
{
    ;
}
void SJRCSBase::proj_and_rebuild(int index,const arma::vec& vox_func,arma::vec& re_vox_func)
{
    ;
}
void SJRCSBase::to_pix_func(int index,const arma::vec& vox_func,arma::vec& pix_func)
{
    ;
}
void SJRCSBase::to_model_func(const arma::mat& alpha,const arma::vec& frame_func,arma::vec& model_func)
{
    ;
}
double SJRCSBase::res_energy(const arma::vec& func0,const arma::vec& func1)
{
    assert(func0.size()==func1.size());
    return arma::accu(arma::square(func0 - func1));
}
void SJRCSBase::max_cir_cor(const arma::vec&,const arma::vec&,arma::uword& max_cir)
{
    ;
}
void SJRCSBase::cir_mat(const arma::sword& cir,arma::mat& alpha)
{
    ;
}
void SJRCSBase::calc_weighted(
        const arma::fmat&vv,
        const arma::fmat&vn,
        arma::Mat<uint8_t>&vc,
        const arma::mat& alpha,
        arma::fmat&wv,
        arma::fmat&wn,
        arma::Mat<uint8_t>&wc
        )
{
    arma::fmat fwc = arma::conv_to<arma::fmat>::from(wc);
    arma::rowvec alpha_median = arma::median( alpha );
    arma::uvec closest_i(alpha.n_cols);

    for(int c=0;c<alpha.n_cols;++c)
    {
        arma::vec col = alpha.col(c);
        arma::uword m;
        col.max(m);
        closest_i(c)=m;
    }

    arma::mat trunc_alpha = alpha;
    for(int c=0;c<alpha.n_cols;++c)
    {
        arma::vec col = trunc_alpha.col(c);
        col( col < alpha_median(c) ).fill(0.0);
        trunc_alpha.col(c) = col;
    }

    arma::rowvec trunc_alpha_colsum = arma::sum(trunc_alpha);

    wv = vv*arma::conv_to<arma::fmat>::from(trunc_alpha);
    wn = vn*arma::conv_to<arma::fmat>::from(trunc_alpha);
    fwc = arma::conv_to<arma::fmat>::from(vc)*arma::conv_to<arma::fmat>::from(trunc_alpha);

    for(int c=0;c<alpha.n_cols;++c)
    {
        if( 0 != trunc_alpha_colsum(c) )
        {
            wv.col(c) /= trunc_alpha_colsum(c);
            wn.col(c) /= trunc_alpha_colsum(c);
            fwc.col(c) /= trunc_alpha_colsum(c);
        }
    }

    wn = arma::normalise( wn );
    wc = arma::conv_to<arma::Mat<uint8_t>>::from(fwc);
}
}
