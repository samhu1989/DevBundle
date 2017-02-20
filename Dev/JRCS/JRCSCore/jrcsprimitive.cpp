#include "jrcsprimitive.h"
#include <QCoreApplication>
#include <QThread>
namespace JRCS{

Plate::Plate():corners_(3,4,arma::fill::zeros),R_(3,3,arma::fill::eye),t_(3,arma::fill::zeros)
{
    xv_.reset(new arma::fmat(3,4));
    xn_.reset(new arma::fmat(3,4));
    xc_.reset(new arma::Mat<uint8_t>(3,4));
    t_ = obj_pos_;
    scale_r_ = arma::linspace<arma::fvec>(0.2,1.2,10);
    trans_r_ = arma::linspace<arma::fvec>(0.2,1.2,10);
}

Plate::Plate(
    const arma::fmat& v,
    const arma::fmat& n,
    const arma::Mat<uint8_t>& c,
    const arma::fvec& pos
):corners_(3,4,arma::fill::zeros),R_(3,3,arma::fill::eye),t_(3,arma::fill::zeros)
{
    xv_.reset(new arma::fmat((float*)v.memptr(),v.n_rows,v.n_cols,false,true));
    xn_.reset(new arma::fmat((float*)n.memptr(),n.n_rows,n.n_cols,false,true));
    xc_.reset(new arma::Mat<uint8_t>((uint8_t*)c.memptr(),c.n_rows,c.n_cols,false,true));
//    std::cerr<<*xv_<<std::endl;
    centroid_ = arma::mean(*xv_,1);
//    std::cerr<<centroid_<<std::endl;
    origin_ = centroid_;
    corners_ = xv_->each_col() - centroid_;
    t_ = arma::fvec(3,arma::fill::zeros);
    size_ = arma::max(arma::abs(corners_),1);
//    std::cerr<<"size_:"<<size_<<std::endl;
    obj_pos_ = pos;
    scale_r_ = arma::linspace<arma::fvec>(0.2,1.2,10);
    trans_r_ = arma::linspace<arma::fvec>(0.2,1.2,10);
}

void Plate::print(void)
{
    std::cerr<<"Plate:"<<std::endl;
    std::cerr<<"xv:"<<std::endl;
    std::cerr<<*xv_<<std::endl;
    std::cerr<<"t:"<<std::endl;
    std::cerr<<t_<<std::endl;
}


void Plate::get_local_translate(
        arma::fvec& t
        )
{
    t = centroid_ - t_;
    t = R_.i()*t;
    t -= origin_;
}

void Plate::local_translate(
        const arma::fvec& t,
        Plate& result
        )
{
//    std::cerr<<"local_translate:  "<<std::endl;
//    std::cerr<<t<<std::endl;
    result.origin_ = origin_;
    *result.xv_ = *xv_;
    //transform back to local coord
    result.xv_->each_col() -= t_;
    *result.xv_ = R_.i()*(*result.xv_);
    //tranlate locally
    result.xv_->each_col() += t;
    //tranlate back
    *result.xv_ = R_*(*result.xv_);
    result.xv_->each_col() += t_;
    result.t_ = R_*t + t_;
    //update centroid
    result.centroid_ = arma::mean(*result.xv_,1);
    result.corners_ = result.xv_->each_col() - result.centroid_;
    if(this!=&result)
    {
        result.R_ = R_;
        *result.xc_ = *xc_;
        *result.xn_ = *xn_;
        result.size_ = size_;
        result.corners_ = corners_;
        result.weighted_centroid_ = weighted_centroid_;
        result.obj_pos_ = obj_pos_;
    }else{
//        std::cerr<<"translate in place"<<std::endl;
    }
}

void Plate::translate(
        const arma::fvec& t,
        Plate& result
        )
{
    result.origin_ = origin_;
    result.t_ = t_ + t;
    *result.xv_ = *xv_;
    result.xv_->each_col() += t;
    //update centroid
    result.centroid_ = arma::mean(*result.xv_,1);
    result.corners_ = result.xv_->each_col() - result.centroid_;
    if(this!=&result)
    {
        result.R_ = R_;
        *result.xc_ = *xc_;
        *result.xn_ = *xn_;
        result.size_ = size_;
        result.weighted_centroid_ = weighted_centroid_;
        result.obj_pos_ = obj_pos_;
    }else{
//        std::cerr<<"translate in place"<<std::endl;
    }
}

void Plate::transform(
        const arma::fmat& R,
        const arma::fvec& t,
        Plate& result
        )
{
    result.origin_ = origin_;
    result.R_ = R*R_;
    result.t_ = R*t_ + t;
    *result.xv_ = R*(*xv_);
    result.xv_->each_col() += t;
    *result.xn_ = R*(*xn_);
    //update centroid
    result.centroid_ = arma::mean(*result.xv_,1);
    result.corners_ = result.xv_->each_col() - result.centroid_;
    if(this!=&result)
    {
        *result.xc_ = *xc_;
        result.size_ = size_;
        result.weighted_centroid_ = weighted_centroid_;
        result.obj_pos_ = obj_pos_;
    }else{
//        std::cerr<<"transform in place"<<std::endl;
    }
}

void Plate::scale(
        const arma::fvec& s,
        Plate& result
        )
{
//    std::cerr<<"scaling:"<<std::endl;
//    std::cerr<<s<<std::endl;
    result.origin_ = origin_;
    result.size_ = size_ % s;
    result.corners_ = R_.i()*corners_;
    result.corners_.each_col() %= s;
    result.corners_ = R_*result.corners_;
    *result.xv_ = result.corners_;
    result.xv_ -> each_col() += centroid_;
    //updating centroid
    result.centroid_ = arma::mean(*result.xv_,1);
    if(this!=&result)
    {
        result.t_ = t_;
        result.R_ = R_;
        *result.xc_ = *xc_;
        *result.xn_ = *xn_;
        result.weighted_centroid_ = weighted_centroid_;
        result.obj_pos_ = obj_pos_;
    }else{
//        std::cerr<<"scale in place"<<std::endl;
    }
}

arma::vec Plate::get_dist2(
        const arma::fmat& v
        )
{
    arma::fmat tv = v.each_col() - t_;
    arma::fmat invR = R_.i();
    assert(invR.is_finite());
    tv = R_.i()*tv;
    return dist(tv,0)+dist(tv,1)+dist(tv,2);
}

arma::vec Plate::dist(const arma::fmat& v, arma::uword dim)
{
    arma::vec result(v.n_cols,arma::fill::zeros);
    arma::frowvec vdim = v.row(dim);
    arma::uvec idx_dim = arma::find( arma::abs( vdim - origin_(dim) ) > size_(dim));
    result(idx_dim) = arma::square( arma::abs( arma::conv_to<arma::vec>::from(vdim.cols(idx_dim)) - origin_(dim) ) - size_(dim) );
    return result;
}

void Plate::get_weighted_centroid(
        const arma::fmat& v,
        const arma::vec& alpha
        )
{
    arma::vec tmp = v * alpha;
    tmp /= arma::accu(alpha);
    weighted_centroid_ = arma::conv_to<arma::fmat>::from(tmp);
}

void Plate::accumulate(
        const arma::fmat& v,
        const arma::fmat& n,
        const arma::Mat<uint8_t>& c,
        const arma::vec alpha
        )
{
    arma::fvec scale_size(3,arma::fill::zeros);
    if(param_.empty())param_ = arma::fcube(scale_r_.size(),scale_r_.size(),trans_r_.size(),arma::fill::zeros);
    else param_.fill(0.0);
    int dim = -1;
    for(int m=0;m<3;++m)
    {
        if(0.0==size_(m))dim=m;
    }
    for(int i = 0;i<scale_r_.size();++i)
    {
        for(int j=0;j<scale_r_.size();++j)
        {
//            #pragma omp parallel for
            for(int k=0;k<trans_r_.size();++k)
            {
                arma::fmat tmpv((float*)xv_->memptr(),xv_->n_rows,xv_->n_cols,true,true);
                arma::fmat tmpn((float*)xn_->memptr(),xn_->n_rows,xn_->n_cols,true,true);
                arma::Mat<uint8_t> tmpc((uint8_t*)xc_->memptr(),xc_->n_rows,xc_->n_cols,true,true);
                Plate::Ptr tmp_plate(new Plate(tmpv,tmpn,tmpc,obj_pos_));
                if(dim==0)
                {
                    scale_size(1) = scale_r_(i);
                    scale_size(2) = scale_r_(j);
                }else if(dim==1)
                {
                    scale_size(0) = scale_r_(i);
                    scale_size(2) = scale_r_(j);
                }else if(dim==2)
                {
                    scale_size(0) = scale_r_(i);
                    scale_size(1) = scale_r_(j);
                }
                scale(scale_size,*tmp_plate);
                arma::fvec t(3,arma::fill::ones);
                arma::fvec t0(3,arma::fill::ones);
                get_local_translate(t);
                t0 = t;
                t(dim) *= trans_r_(k);
//                std::cerr<<"t in accumulate"<<t<<std::endl;
                tmp_plate->local_translate((t-t0),*tmp_plate);
                arma::vec dist2 = tmp_plate->get_dist2(v);
                dist2 %= alpha;
                param_(i,j,k) = arma::accu(dist2);
            }
        }
        QCoreApplication::processEvents();
    }
}

void Plate::accumulate(const Plate& p)
{
    if(param_.empty())
    {
        param_ = p.param_;
    }else{
        param_ += p.param_;
    }
}

void Plate::fit(void)
{
    arma::uword i,j,k;
    param_.min(i,j,k);
    //find the minimum
//    std::cerr<<"the minimum:"<<i<<","<<j<<","<<k<<std::endl;
    //update this with the minimum
    int dim = -1;
    for(int m=0;m<3;++m)
    {
        if(0.0==size_(m))dim=m;
    }
    assert(dim>=0&&dim<=2);
    arma::fvec scale_size(3,arma::fill::zeros);
    if(dim==0)
    {
        scale_size(1) = scale_r_(i);
        scale_size(2) = scale_r_(j);
    }else if(dim==1)
    {
        scale_size(0) = scale_r_(i);
        scale_size(2) = scale_r_(j);
    }else if(dim==2)
    {
        scale_size(0) = scale_r_(i);
        scale_size(1) = scale_r_(j);
    }
    scale(scale_size,*this);
    arma::fvec t(3,arma::fill::ones);
    arma::fvec t0(3,arma::fill::ones);
    get_local_translate(t);
    t0 = t;
    t(dim) *= trans_r_(k);
    local_translate( (t - t0) , *this );
    param_.clear();
}

JRCSPrimitive::JRCSPrimitive():JRCSBilateral(),plate_num_for_obj_(5),point_num_for_plate_(4)
{

}

void JRCSPrimitive::initx(
        const MatPtr& xv,
        const MatPtr& xn,
        const CMatPtr& xc
        )
{
    if(verbose_)std::cerr<<"SJRCSBase start init x"<<std::endl;

    xv_ptr_ = xv;
    xn_ptr_ = xn;
    xc_ptr_ = xc;

    max_obj_radius_ = 5.0;// for randomly reset rt

    xv_ptr_->fill(std::numeric_limits<float>::quiet_NaN());
    xn_ptr_->fill(std::numeric_limits<float>::quiet_NaN());
    xc_ptr_->fill(0);

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
    plate_ptrlst_.resize(plate_num_for_obj_*obj_num_);
    arma::uword* pxr = (arma::uword*)obj_range_.data();
    arma::uword* pxr_s = pxr;
    *pxr = 0;
    for(int obj_idx = 0 ; obj_idx < obj_num_ ; ++ obj_idx )
    {
        int obj_size = plate_num_for_obj_*point_num_for_plate_;
        if( pxr > pxr_s )*pxr = (*(pxr - 1))+1;
        *(pxr+1) = *pxr + plate_num_for_obj_ - 1;

        arma::fmat objv(pxv,3,obj_size,false,true);
        arma::fmat objn(pxn,3,obj_size,false,true);
        arma::Mat<uint8_t> objc(pxc,3,obj_size,false,true);

        arma::fvec pos = obj_pos_.col(obj_idx);
        reset_obj_vn(0.05,pos,objv,objn);
        reset_obj_c(objc);

        for(int j=0 ; j < plate_num_for_obj_; ++j)
        {
            plate_ptrlst_[obj_idx*plate_num_for_obj_+j].reset(
                new Plate(
                    arma::fmat(pxv,3,point_num_for_plate_,false,true),
                    arma::fmat(pxn,3,point_num_for_plate_,false,true),
                    arma::Mat<uint8_t>(pxc,3,point_num_for_plate_,false,true),
                    obj_pos_.col(obj_idx)
                )
            );
            pxv += 3*point_num_for_plate_;
            pxn += 3*point_num_for_plate_;
            pxc += 3*point_num_for_plate_;
        }

        pxr += 2;
        r_k -= obj_size;
    }
    if(verbose_)std::cerr<<"Done initx_"<<std::endl;
}

void JRCSPrimitive::reset_obj_vn(
        float radius,
        arma::fvec& pos,
        arma::fmat& ov,
        arma::fmat& on
        )
{
    arma::fmat v = {
        {-1, 1, 1,-1, 1,-1,-1, 1, 1, 1, 1, 1,-1,-1,-1,-1,-1,-1, 1, 1},
        { 1, 1, 1, 1,-1,-1,-1,-1, 1,-1,-1, 1,-1, 1, 1,-1, 1,-1,-1, 1},
        { 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 1, 1}
    };
    ov = v;
    ov *= radius;
    ov.each_col() += pos;
    arma::fmat n = {
        { 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1,-1,-1,-1,-1, 0, 0, 0, 0},
        { 1, 1, 1, 1,-1,-1,-1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1}
    };
    on = n;
}

void JRCSPrimitive::compute(void)
{
    if(verbose_)std::cerr<<"preparing"<<std::endl;
    prepare_primitive();
    while(!isEnd_primitive())
    {
        update_color_label();
        if(verbose_)std::cerr<<"step 1"<<std::endl;
        #pragma omp parallel for
        for( int i=0 ; i < vvs_ptrlst_.size() ; ++i )
        {
            step_1(i);
        }
        if(verbose_)std::cerr<<"step 2"<<std::endl;
        step_2();
        finish_primitive();
//        QThread::currentThread()->sleep(60);
    }
//    JRCSBilateral::compute();
}

void JRCSPrimitive::prepare_primitive()
{
    iter_count_ = 0;
    reset_alpha_primitive();
    plate_t_ptrlst_.resize(vvs_ptrlst_.size());
    #pragma omp parallel for
    for( int i = 0 ; i < vvs_ptrlst_.size() ; ++i )
    {
        plate_t_ptrlst_[i].resize(plate_ptrlst_.size());
        float* pwv = (float*)wvs_ptrlst_[i]->memptr();
        float* pwn = (float*)wns_ptrlst_[i]->memptr();
        uint8_t* pwc = (uint8_t*)wcs_ptrlst_[i]->memptr();
        for(int j=0 ; j < plate_t_ptrlst_[i].size() ; ++j )
        {
            plate_t_ptrlst_[i][j].reset(
                        new Plate(
                            arma::fmat(pwv,3,point_num_for_plate_,false,true),
                            arma::fmat(pwn,3,point_num_for_plate_,false,true),
                            arma::Mat<uint8_t>(pwc,3,point_num_for_plate_,false,true),
                            arma::fvec(3,arma::fill::zeros)
                            )
                        );
            pwv += 3*point_num_for_plate_;
            pwn += 3*point_num_for_plate_;
            pwc += 3*point_num_for_plate_;
        }
    }
    reset_prob_primitive();
}

void JRCSPrimitive::finish_primitive()
{
    ++iter_count_;
}

void JRCSPrimitive::step_1(int i)
{
    //input
    arma::fmat& vv_ = *vvs_ptrlst_[i];
    arma::fmat& vn_ = *vns_ptrlst_[i];
    arma::Mat<uint8_t>& vc_ = *vcs_ptrlst_[i];
    if(verbose_>1)std::cerr<<"transform latent model"<<std::endl;
    Ts& rt = rt_lst_[i];
    for(int o = 0 ; o < obj_num_ ; ++o )
    {
        arma::fmat R(rt[o].R,3,3,false,true);
        arma::fvec t(rt[o].t,3,false,true);

        for(int j=obj_range_[2*o];j<=obj_range_[2*o+1];++j)
        {
            plate_ptrlst_[j]->transform(R,t,*plate_t_ptrlst_[i][j]);
        }
    }

    if(verbose_>1)std::cerr<<"calculate alpha"<<std::endl;
    arma::mat&  alpha = *alpha_ptrlst_[i];
    for(int c = 0 ; c < alpha.n_cols ; ++c )
    {
        alpha.col(c) = plate_t_ptrlst_[i][c]->get_dist2(vv_);
    }

    alpha.each_row() %= (-0.5*xv_invvar_);
    alpha = arma::trunc_exp(alpha);

    alpha.each_row() %=  arma::pow(xv_invvar_,1.5);

    if(verbose_>1)std::cerr<<"normalize alpha"<<std::endl;
    alpha += std::numeric_limits<double>::epsilon(); //add eps for numeric stability
    arma::vec alpha_rowsum = ( 1.0 + beta_ ) * arma::sum(alpha,1);
    alpha.each_col() /= alpha_rowsum;

    // cut the ones below median to zeros for better converge
    for(int c = 0; c < alpha.n_cols ; ++c )
    {
        arma::vec col = alpha.col(c);
        col( col < arma::median(col) ).fill(0.0);
        alpha.col(c) = col;
    }

    if(verbose_>1)std::cerr<<"normalize alpha"<<std::endl;
    alpha += std::numeric_limits<double>::epsilon(); //add eps for numeric stability
    alpha_rowsum = ( 1.0 + beta_ ) * arma::sum(alpha,1);
    alpha.each_col() /= alpha_rowsum;

    if(!alpha.is_finite())
    {
        std::cerr<<iter_count_<<":invalid alpha in step_a("<<i<<")"<<std::endl;
    }

    //update RT
    if(verbose_>1)std::cerr<<"#1 calculate weighted plate centers"<<std::endl;
    arma::frowvec alpha_colsum = arma::conv_to<arma::frowvec>::from(arma::sum( alpha ));
    assert(plate_t_ptrlst_[i].size()==alpha.n_cols);
    for(int c = 0 ; c < alpha.n_cols ; ++c )
    {
        plate_t_ptrlst_[i][c]->get_weighted_centroid(vv_,alpha.col(c));
//        if(i==0&&0==c%5)
//        {
//            std::cerr<<"c:"<<plate_t_ptrlst_[i][c]->centroid_.t()<<std::endl;
//            std::cerr<<"wc:"<<plate_t_ptrlst_[i][c]->weighted_centroid_.t()<<std::endl;
//        }
    }

    if(verbose_>1)std::cerr<<"#2 calculating RT for each object"<<std::endl;
    for(int o = 0 ; o < obj_num_ ; ++o )
    {
        arma::fmat A;
        arma::fmat U,V;
        arma::fvec s;
        arma::fmat dR;
        arma::fvec dt;
        arma::fmat R(rt[o].R,3,3,false,true);
        arma::fvec t(rt[o].t,3,false,true);
        arma::fmat _v(3,plate_num_for_obj_);
        arma::fmat objv(3,plate_num_for_obj_);
        int vi = 0;
        for(int p=obj_range_[2*o];p<=obj_range_[2*o+1];++p)
        {
            _v.col(vi) = plate_t_ptrlst_[i][p]->weighted_centroid_;
            objv.col(vi) = plate_t_ptrlst_[i][p]->centroid_;
            ++vi;
        }
        arma::fmat cv = _v.each_col() - arma::mean(_v,1);
        if(!_v.is_finite())
        {
            std::cerr<<iter_count_<<":"<<o<<":!v.is_finite()"<<std::endl;
        }

        arma::frowvec square_lambdav = arma::conv_to<arma::frowvec>::from(xv_invvar_.cols(obj_range_[2*o],obj_range_[2*o+1]));
        square_lambdav %= alpha_colsum.cols(obj_range_[2*o],obj_range_[2*o+1]);
        square_lambdav += std::numeric_limits<float>::epsilon(); // for numeric stability
        if(!square_lambdav.is_finite())
        {
            std::cerr<<"!square_lambda.is_finite()"<<std::endl;
        }

        arma::frowvec pv(square_lambdav.n_cols,arma::fill::ones);
        arma::frowvec square_norm_lambdav = square_lambdav / arma::accu(square_lambdav);
        if(!square_norm_lambdav.is_finite())
        {
            std::cerr<<"!square_norm_lambda.is_finite()"<<std::endl;
        }
        pv -= square_norm_lambdav;
        arma::fmat tmpv = objv.each_col() - arma::mean(objv,1) ;
        tmpv.each_row() %= pv % square_lambdav;

        A = cv*tmpv.t();

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
                arma::fmat tmp = _v - dR*objv;
                tmp.each_row() %= square_norm_lambdav;
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
                arma::fmat tmp = _v - dR*objv;
                tmp.each_row() %= square_norm_lambdav;
                dt = arma::sum(tmp,1);
                if(!dt.is_finite())
                {
                    std::cerr<<iter_count_<<":!dt.is_finite()"<<std::endl;
                    dt.fill(arma::fill::zeros);
                }
            }
        }
        }
        //transforming transformed object
        for(int p=obj_range_[2*o];p<=obj_range_[2*o+1];++p)
        {
            plate_t_ptrlst_[i][p]->transform(dR,dt,*plate_t_ptrlst_[i][p]);
        }

        //updating R T
        R = dR*R;
        t = dR*t + dt;
    }
    if(verbose_>1)std::cerr<<"#3 done RT for each object"<<std::endl;
    arma::mat alpha_v2(alpha.n_rows,alpha.n_cols);
    for(int c=0;c<alpha_v2.n_cols;++c)
    {
        alpha_v2.col(c) = plate_t_ptrlst_[i][c]->get_dist2(vv_);
    }
    vvar_.row(i) = arma::sum(alpha_v2%alpha);

    if(verbose_>1)std::cerr<<"#4 done var for each object"<<std::endl;
    for(int c=0;c<alpha.n_cols;++c)
    {
        plate_t_ptrlst_[i][c]->accumulate(vv_,vn_,vc_,alpha.col(c));
    }
    if(verbose_>1)std::cerr<<"#5 done accumulation"<<std::endl;
}

void JRCSPrimitive::step_2(void)
{
    if(verbose_>1)std::cerr<<"Updating Latent Model"<<std::endl;
    #pragma omp parallel for
    for( int i=0 ; i < plate_ptrlst_.size() ; ++i )
    {
        for(int j=0;j<plate_t_ptrlst_.size();++j)
        {
            plate_ptrlst_[i]->accumulate(*plate_t_ptrlst_[j][i]);
        }
        plate_ptrlst_[i]->fit();
    }

    if(verbose_>1)std::cerr<<"Updating variance of Latent Model"<<std::endl;
    arma::rowvec alpha_sum;
    for(int i=0;i<alpha_ptrlst_.size();++i)
    {
        arma::mat& alpha = *alpha_ptrlst_[i];
        if(alpha_sum.empty())alpha_sum = arma::sum(alpha);
        else alpha_sum += arma::sum(alpha);
    }

    if(verbose_>1)std::cerr<<"Restore reciprocal fractions of variation"<<std::endl;
    xv_invvar_ = ( ( 3.0*alpha_sum ) / ( arma::sum(vvar_) + beta_ ) );

    if(verbose_>1)std::cerr<<"Updating weight of centroids"<<std::endl;
    float mu = arma::accu(alpha_sum);
    mu *= ( 1.0 + beta_ );
    x_p_ = alpha_sum;
    if( mu != 0)x_p_ /= mu;
}

bool JRCSPrimitive::isEnd_primitive(void)
{
    if(iter_count_>=max_init_iter_)return true;
    else return false;
}

void JRCSPrimitive::reset_alpha_primitive()
{
    if(verbose_>0)std::cerr<<"allocating alpha"<<std::endl;
    int N = plate_ptrlst_.size();
    int idx=0;
    while( idx < vvs_ptrlst_.size() )
    {
        if(idx>=alpha_ptrlst_.size())alpha_ptrlst_.emplace_back(new arma::mat(vvs_ptrlst_[idx]->n_cols,N));
        else if((alpha_ptrlst_[idx]->n_rows!=vvs_ptrlst_[idx]->n_cols)||(alpha_ptrlst_[idx]->n_cols!=N))
            alpha_ptrlst_[idx].reset(new arma::mat(vvs_ptrlst_[idx]->n_cols,N));
        ++idx;
    }
    if(verbose_>0)std::cerr<<"done allocating alpha"<<std::endl;
}

void JRCSPrimitive::reset_prob_primitive()
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
    xv_invvar_ = arma::rowvec(plate_ptrlst_.size());
    xv_invvar_.fill(1.0/maxvar);

    vvar_ = arma::mat(vvs_ptrlst_.size(),plate_ptrlst_.size());

    x_p_ = arma::rowvec(plate_ptrlst_.size());
    x_p_.fill(1.0/float(plate_ptrlst_.size()));

    if(verbose_>0)std::cerr<<"done probability"<<std::endl;
}

}
