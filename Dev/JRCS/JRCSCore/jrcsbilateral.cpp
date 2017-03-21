#include "jrcsbilateral.h"
namespace JRCS{
JRCSBilateral::JRCSBilateral():SJRCSBase()
{
    ;
}
//Step A:
//Calculate Alpha(Expectation)
//Update RT
//Step B:
//Update X
//Update var pk
void JRCSBilateral::compute(void)
{

    prepare_compute();
    while(!isEnd())
    {
        if(verbose_)std::cerr<<"step a"<<std::endl;
        #pragma omp parallel for
        for( int i=0 ; i < vvs_ptrlst_.size() ; ++i )
        {
            step_a(i);
        }
        if(verbose_)std::cerr<<"step b"<<std::endl;
        step_b();
        finish_steps();
    }
}

void JRCSBilateral::prepare_compute(void)
{
    if(verbose_)std::cerr<<"JRCSBilateral::preparing"<<std::endl;
    iter_count_ = 0;
    reset_alpha();
    xtv_ptrlst_.resize(vvs_ptrlst_.size());
    xtn_ptrlst_.resize(vns_ptrlst_.size());
    #pragma omp for
    for( int i = 0 ; i < vvs_ptrlst_.size() ; ++i )
    {
        xtv_ptrlst_[i].reset(new arma::fmat(xv_ptr_->n_rows,xv_ptr_->n_cols,arma::fill::zeros));
        xtn_ptrlst_[i].reset(new arma::fmat(xn_ptr_->n_rows,xn_ptr_->n_cols,arma::fill::zeros));
    }
    reset_prob();
    if(!xf_invvar_.is_finite())
    {
        std::cerr<<"invalid xf_invvar_"<<std::endl;
    }
    if(!xf_ptr_->is_finite())
    {
        std::cerr<<"invalid xf_ptr_ in prepare"<<std::endl;
    }
    rescale_feature();
    QCoreApplication::processEvents();
}

void JRCSBilateral::step_a(int i)
{
    //latent model
    arma::fmat& xv_ = *xv_ptr_;
    arma::fmat& xn_ = *xn_ptr_;
    arma::fmat& xf_ = *xf_ptr_;

    //input
    arma::fmat& vv_ = *vvs_ptrlst_[i];
    arma::fmat& vn_ = *vns_ptrlst_[i];
    arma::fmat& vf_ = *vfs_ptrlst_[i];
    if(verbose_>1)std::cerr<<"transform latent model"<<std::endl;
    arma::fmat& xtv_ = *xtv_ptrlst_[i];
    arma::fmat& xtn_ = *xtn_ptrlst_[i];
    Ts& rt = rt_lst_[i];
    for(int o = 0 ; o < obj_num_ ; ++o )
    {
        arma::fmat R(rt[o].R,3,3,false,true);
        arma::fvec t(rt[o].t,3,false,true);
        arma::uword offset = 3*obj_range_[2*o];
        arma::uword size = obj_range_[2*o+1] - obj_range_[2*o] + 1;
        arma::fmat objv(((float*)xtv_.memptr())+offset,3,size,false,true);
        arma::fmat objn(((float*)xtn_.memptr())+offset,3,size,false,true);
        objv = R*xv_.cols(obj_range_[2*o],obj_range_[2*o+1]);
        objv.each_col() += t;
        objn = R*xn_.cols(obj_range_[2*o],obj_range_[2*o+1]);
    }
    if(verbose_>1)std::cerr<<"calculate alpha"<<std::endl;

    arma::mat&  alpha = *alpha_ptrlst_[i];
    double dim = double(f_dim_)/2.0f;
    for(int r = 0 ; r < alpha.n_rows ; ++r )
    {
        arma::mat tmpf = arma::conv_to<arma::mat>::from( xf_.each_col() - vf_.col(r) );
        arma::mat tmpv = arma::conv_to<arma::mat>::from( xtv_.each_col() - vv_.col(r) );

        arma::rowvec alpha_f = arma::sum(arma::square(tmpf))%(-0.5*xf_invvar_) ;
        arma::rowvec alpha_v = arma::sum(arma::square(tmpv))%(-0.5*xv_invvar_);

        alpha_f = arma::trunc_exp(alpha_f);
        alpha_f %= arma::pow(xf_invvar_,dim);

        alpha_v = arma::trunc_exp(alpha_v);
        alpha_v %= arma::pow(xv_invvar_,1.5);

        alpha.row(r) = alpha_v % alpha_f ;
    }

    if(verbose_>1)std::cerr<<"normalize alpha"<<std::endl;
    alpha += std::numeric_limits<double>::epsilon(); //add eps for numeric stability
    arma::vec alpha_rowsum = ( 1.0 + beta_ ) * arma::sum(alpha,1);
    alpha.each_col() /= alpha_rowsum;

    if(!alpha.is_finite())
    {
        std::cerr<<iter_count_<<":invalid alpha in step_a("<<i<<")"<<std::endl;
    }

    //update RT
    if(verbose_>1)std::cerr<<"#1 calculate weighted point cloud"<<std::endl;
    arma::Mat<uint8_t>& vc = *vcs_ptrlst_[i];
    arma::fmat& wv = *wvs_ptrlst_[i];
    arma::fmat& wn = *wns_ptrlst_[i];
    arma::Mat<uint8_t>& wc = *wcs_ptrlst_[i];
    arma::fmat& wf = *wfs_ptrlst_[i];
    arma::frowvec alpha_colsum = arma::conv_to<arma::frowvec>::from(arma::sum( alpha ));
    calc_weighted(vv_,vn_,vf_,vc,alpha,wv,wn,wf,wc);

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
        arma::uword offset = 3*obj_range_[2*o];
        arma::uword size = obj_range_[2*o+1] - obj_range_[2*o] + 1;
        arma::fmat objv(((float*)xtv_.memptr())+offset,3,size,false,true);
        arma::fmat objn(((float*)xtn_.memptr())+offset,3,size,false,true);
        arma::fmat _v = wv.cols(obj_range_[2*o],obj_range_[2*o+1]);
        arma::fmat _n = wn.cols(obj_range_[2*o],obj_range_[2*o+1]);
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

//        std::cerr<<"od:"<<o<<std::endl;

        //transforming transformed object
        objv = dR*objv;
        objv.each_col() += dt;

        objn = dR*objn;

        //updating R T
        R = dR*R;
        t = dR*t + dt;
    }
    if(verbose_>1)std::cerr<<"#3 done RT for each object"<<std::endl;

    arma::fmat alpha_v2(alpha.n_rows,alpha.n_cols);
    for(int r=0;r<alpha_v2.n_rows;++r)
    {
        alpha_v2.row(r) = arma::sum(arma::square(xtv_.each_col() - vv_.col(r)));
    }
    vvar_.row(i) = arma::sum(alpha_v2%alpha);

    arma::fmat alpha_f2(alpha.n_rows,alpha.n_cols);
    for(int r=0;r<alpha_f2.n_rows;++r)
    {
        alpha_f2.row(r) = arma::sum(arma::square(xf_.each_col() - vf_.col(r)));
    }
    fvar_.row(i) = arma::sum(alpha_f2%alpha);
}

void JRCSBilateral::step_b(void)
{
    xv_ptr_->fill(0.0);
    xn_ptr_->fill(0.0);
    xf_ptr_->fill(0.0);
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
        arma::fmat& wf = *wfs_ptrlst_[i];
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
            xf_ptr_->cols(s,e) = rate0*xf_ptr_->cols(s,e)+rate1*wf.cols(s,e);
            xc.cols(s,e) = rate0*xc.cols(s,e)+rate1*wc.cols(s,e);
        }
    }

    //the method of updating normal is different from others
    int o = 0;
    for( int i=0 ; i < xn_ptr_->n_cols ; ++i )
    {
        if( i > obj_range_[2*o+1] ) ++o;
        arma::mat X(3,wns_ptrlst_.size());
        #pragma omp parallel for
        for( int f = 0 ; f < wns_ptrlst_.size() ; ++ f )
        {
            Ts& rt = rt_lst_[f];
            arma::fmat R(rt[o].R,3,3,false,true);
            arma::fmat invR = R.i();
            X.col(f) = arma::conv_to<arma::vec>::from( invR*wns_ptrlst_[f]->col(i) );
        }
        arma::mat eigvec;
        arma::vec eigval;
        arma::eig_sym(eigval,eigvec,(X*X.t()),"std");
        xn_ptr_->col(i) = arma::conv_to<arma::fvec>::from(eigvec.col(2));
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

    if(verbose_>1)std::cerr<<"Fix the x center position"<<std::endl;
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

    if(verbose_>1)std::cerr<<"Updating variance of centroids"<<std::endl;
    arma::rowvec alpha_sum;
    for(int i=0;i<alpha_ptrlst_.size();++i)
    {
        arma::mat& alpha = *alpha_ptrlst_[i];
        if(alpha_sum.empty())alpha_sum = arma::sum(alpha);
        else alpha_sum += arma::sum(alpha);
    }

    if(verbose_>1)std::cerr<<"Restore reciprocal fractions of variation"<<std::endl;
    xv_invvar_ = ( ( 3.0*alpha_sum ) / ( arma::sum(vvar_) + beta_ ) );
    xf_invvar_ = ( ( double(f_dim_)*alpha_sum ) / ( arma::sum(fvar_) + beta_ ) );

    if(verbose_>1)std::cerr<<"Updating weight of centroids"<<std::endl;
    float mu = arma::accu(alpha_sum);
    mu *= ( 1.0 + beta_ );
    x_p_ = alpha_sum;
    if( mu != 0)x_p_ /= mu;
}

void JRCSBilateral::reset_prob()
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
    xv_invvar_ = arma::rowvec(xv_ptr_->n_cols);
    xv_invvar_.fill(1.0/maxvar);

    xf_invvar_ = arma::rowvec(xv_ptr_->n_cols);
    xf_invvar_.fill(2.0);

    vvar_ = arma::mat(vvs_ptrlst_.size(),xv_ptr_->n_cols);
    fvar_ = arma::mat(vvs_ptrlst_.size(),xv_ptr_->n_cols);

    x_p_ = arma::rowvec(xv_ptr_->n_cols);
    x_p_.fill(1.0/float(xv_ptr_->n_cols));

    if(verbose_>0)std::cerr<<"done probability"<<std::endl;
}

void JRCSBilateral::calc_obj()
{
    obj_v = 0.0;
    obj_f = 0.0;
    for(int idx=0;idx<vvs_ptrlst_.size();++idx)
    {
        arma::mat& alpha = *alpha_ptrlst_[idx];
        arma::fmat& vv_ = *vvs_ptrlst_[idx];
        arma::fmat& vn_ = *vns_ptrlst_[idx];
        arma::fmat& xtv_ = *xtv_ptrlst_[idx];
        arma::fmat& xtn_ = *xtn_ptrlst_[idx];
        arma::fmat alpha_v2(alpha.n_rows,alpha.n_cols);
        #pragma omp parallel for
        for(int r=0;r<alpha_v2.n_rows;++r)
        {
            arma::frowvec tmpv = arma::sum(arma::square(xtv_.each_col() - vv_.col(r)));
            tmpv %= arma::conv_to<arma::frowvec>::from(xv_invvar_);
            tmpv -= 1.5*arma::conv_to<arma::frowvec>::from(arma::trunc_log(xv_invvar_));
            tmpv -= 1.0*arma::conv_to<arma::frowvec>::from(arma::trunc_log(x_p_));
            alpha_v2.row(r) = tmpv;
        }
        obj_v += arma::accu( alpha % alpha_v2 );
        arma::fmat alpha_n2(alpha.n_rows,alpha.n_cols);
        #pragma omp parallel for
        for(int r=0;r<alpha_n2.n_rows;++r)
        {
            arma::frowvec tmpn = arma::sum(arma::square(xtn_.each_col() - vn_.col(r)));
            tmpn %= arma::conv_to<arma::frowvec>::from(xf_invvar_);
            tmpn -= double(f_dim_)/2.0*arma::conv_to<arma::frowvec>::from(arma::trunc_log(xf_invvar_));
            tmpn -= 1.0*arma::conv_to<arma::frowvec>::from(arma::trunc_log(x_p_));
            alpha_n2.row(r) = tmpn;
        }
        obj_f += arma::accu( alpha % alpha_n2 );
    }
    std::cout<<"obj_v:"<<obj_v<<std::endl;
    std::cout<<"obj_f:"<<obj_f<<std::endl;
    obj_vec_.push_back(obj_f + obj_v) ;
}

void JRCSBilateral::calc_weighted(
        const arma::fmat&vv,
        const arma::fmat&vn,
        const arma::fmat&vf,
        arma::Mat<uint8_t>&vc,
        const arma::mat& alpha,
        arma::fmat&wv,
        arma::fmat&wn,
        arma::fmat&wf,
        arma::Mat<uint8_t>&wc
        )
{
    if(verbose_>1)std::cerr<<"JRCSBilateral::calc_weighted"<<std::endl;
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
    wf = vf*arma::conv_to<arma::fmat>::from(trunc_alpha);

    //normal is weighted differently
    for(int c=0;c<wn.n_cols;++c)
    {
        arma::mat wvn = arma::conv_to<arma::mat>::from(vn);
        wvn.each_row() %= trunc_alpha.col(c).t();
        arma::vec eigval;
        arma::mat eigvec;
        arma::mat X = wvn*wvn.t();
        arma::eig_sym(eigval,eigvec,X,"std");
        if( eigvec.n_cols <3 )
        {
            std::cerr<<"eigvec.n_cols:"<<eigvec.n_cols<<std::endl;
            std::cerr<<"X:"<<X<<std::endl;
        }
        wn.col(c) = arma::conv_to<arma::fvec>::from(eigvec.col(2));
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
            wf.col(c) /= trunc_alpha_colsum(c);
            fwc.col(c) /= trunc_alpha_colsum(c);
        }
    }

    if(!wv.is_finite())
    {
        std::cerr<<iter_count_<<":b:!wv.is_finite()"<<std::endl;
    }

    if(!wf.is_finite())
    {
        std::cerr<<iter_count_<<":b:!wf.is_finite()"<<std::endl;
    }

//    wn = arma::normalise( wn );
    assert(wv.is_finite());
    assert(wn.is_finite());
    assert(wf.is_finite());
    assert(fwc.is_finite());
    wc = arma::conv_to<arma::Mat<uint8_t>>::from(fwc);
    if(verbose_>1)std::cerr<<"JRCSBilateral::calc_weighted finished"<<std::endl;
}

bool JRCSBilateral::configure(Config::Ptr config)
{
    if(!SJRCSBase::configure(config))return false;
    use_res_ = false;
    return true;
}

bool JRCSBilateral::input_extra(const MeshBundle<DefaultMesh>::PtrList& inputs)
{
    bool use_feature_ = true;
    f_dim_ = 3;
    foreach(MeshBundle<DefaultMesh>::Ptr m_ptr,inputs)
    {
        if( m_ptr->p_feature_.empty() || m_ptr->p_feature_.n_cols != m_ptr->mesh_.n_vertices())
        {
            use_feature_ = false;
            std::cerr<<"using point color"<<std::endl;
            break;
        }else{
            if( m_ptr->p_feature_.n_rows > f_dim_ )
            {
                f_dim_ = m_ptr->p_feature_.n_rows;
            }
        }
    }
    if(use_feature_){
        std::cerr<<"using external feature"<<std::endl;
    }
    xf_ptr_.reset(new arma::fmat(f_dim_,xv_ptr_->n_cols,arma::fill::zeros));
    vfs_ptrlst_.reserve(inputs.size());
    wfs_ptrlst_.reserve(inputs.size());
    arma::uword i = 0;
    double max_fvar = std::numeric_limits<float>::lowest();
    foreach(MeshBundle<DefaultMesh>::Ptr m_ptr,inputs)
    {
        if(use_feature_)
        {

            vfs_ptrlst_.emplace_back(new arma::fmat(
                                         m_ptr->p_feature_.memptr(),
                                         m_ptr->p_feature_.n_rows,
                                         m_ptr->p_feature_.n_cols,
                                         false,
                                         true
                                         )
                                     );
        }else{

            vfs_ptrlst_.emplace_back(new arma::fmat(3,m_ptr->mesh_.n_vertices()));
            ( *vfs_ptrlst_.back() ) = arma::conv_to<arma::fmat>::from(*vcs_ptrlst_[i]);

        }
        arma::fvec fmean;
        arma::fvec fmin,fmax;
        arma::fmat& fs = *vfs_ptrlst_.back();
        fmin = arma::min(fs,1);
        fmax = arma::max(fs,1);
        fmean = arma::mean(fs,1);
        if(!fmean.is_finite())
        {
            std::cerr<<"invalid fmean"<<std::endl;
        }
        xf_ptr_->each_col() += fmean;
        float fvar = arma::max(arma::square(fmax - fmin));
        if( fvar > max_fvar ) max_fvar = fvar;
        wfs_ptrlst_.emplace_back(new arma::fmat(xf_ptr_->n_rows,xf_ptr_->n_cols));
        ++i;
    }
    (*xf_ptr_) /= float(inputs.size());
    if(!xf_ptr_->is_finite())
    {
        std::cerr<<"invalid xf_ptr_ in input extra"<<std::endl;
    }
    if(verbose_)std::cerr<<"max_fvar:"<<max_fvar<<std::endl;
    xf_invvar_.fill(1.0/max_fvar);
    return true;
}

void JRCSBilateral::rescale_feature()
{
    arma::mat X(f_dim_,f_dim_,arma::fill::zeros);
    arma::fvec fmean(f_dim_);
    for(int i = 0 ; i < vfs_ptrlst_.size() ; i ++ )
    {
        arma::fmat& f = *vfs_ptrlst_[i];
        fmean += arma::mean(f,1);
    }
    float k = vfs_ptrlst_.size();
    fmean /= k;
    double N = 0.0;
    for(int i = 0 ; i < vfs_ptrlst_.size() ; i ++ )
    {
        arma::fmat& f = *vfs_ptrlst_[i];
        f.each_col() -= fmean;
        X += arma::conv_to<arma::mat>::from( f*f.t() );
        N += double(f.n_cols);
    }
    arma::mat eigvec;
    arma::vec eigval;
    arma::eig_sym(eigval,eigvec,X);
    std::cerr<<"eigval:"<<eigval<<std::endl;
    arma::fvec scale(eigval.size());
    scale = arma::conv_to<arma::fvec>::from( 1.0 / arma::sqrt(arma::abs( eigval / N )));
    //rotate feature space
    #pragma omp parallel for
    for(int i = 0 ; i < vfs_ptrlst_.size() ; i ++ )
    {
        arma::fmat& f = *vfs_ptrlst_[i];
        f = ( f.t() * arma::conv_to<arma::fmat>::from(eigvec) ).t();
    }
    xf_invvar_.fill(1.0);
}
}
