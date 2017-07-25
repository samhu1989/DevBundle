#include "jrcsbase.h"
#include <QThread>
#include <strstream>
#include "MeshColor.h"
#include "densecrf3d.h"

namespace JRCS{

bool JRCSBase::configure(Config::Ptr config)
{
    std::cerr<<"configuring"<<std::endl;
    config_ = config;
    if(config_->has("JRCS_obj_w"))
    {
        std::vector<float> objw;
        config_->getFloatVec("JRCS_obj_w",objw);
        reset_objw(objw);
        if(0==obj_num_)return false;
    }else return false;

    if(config_->has("JRCS_max_iter"))
    {
        set_max_iter(config_->getInt("JRCS_max_iter"));
        obj_vec_.reserve(config_->getInt("JRCS_max_iter"));
    }else return false;

    if(config_->has("JRCS_max_init"))
    {
        set_max_init_iter(config_->getInt("JRCS_max_init"));
    }else set_max_init_iter(config_->getInt("JRCS_max_iter")/2);

    if(config_->has("JRCS_verbose"))
    {
        verbose_=config_->getInt("JRCS_verbose");
    }else{
        verbose_=-1;
    }

    if(config_->has("JRCS_smooth"))
    {
        if(!config_->getInt("JRCS_smooth"))enable_smooth(false);
        else{
            enable_smooth(true);
            if(config_->has("JRCS_smooth_w"))
            {
                set_smooth_weight(config_->getFloat("JRCS_smooth_w"));
            }else set_smooth_weight(1.0);
            if(config_->has("JRCS_smooth_iter"))
            {
                set_max_smooth_iter(config_->getInt("JRCS_smooth_iter"));
            }else set_max_smooth_iter(1);
        }
    }else{
        enable_smooth(true);
        if(config_->has("JRCS_smooth_w"))
        {
            set_smooth_weight(config_->getFloat("JRCS_smooth_w"));
        }else set_smooth_weight(1.0);
        if(config_->has("JRCS_smooth_iter"))
        {
            set_max_smooth_iter(config_->getInt("JRCS_smooth_iter"));
        }else set_max_smooth_iter(1);
    }

    if(config_->has("JRCS_debug_path"))
    {
        set_debug_path(config_->getString("JRCS_debug_path"));
    }else set_debug_path("./debug/");

    if(config_->has("JRCS_mu_type"))
    {
        if(config_->getString("JRCS_mu_type")=="ObjOnly")set_mu_type(JRCS::JRCSBase::ObjOnly);
        if(config_->getString("JRCS_mu_type")=="ObjPointDist")set_mu_type(JRCS::JRCSBase::ObjPointDist);
    }else set_mu_type(JRCS::JRCSBase::ObjOnly);

    if(config_->has("JRCS_rt_type"))
    {
        if(config_->getString("JRCS_rt_type")=="Gamma")set_rt_type(JRCS::JRCSBase::Gamma);
    }else set_rt_type(JRCS::JRCSBase::All);

    if(config_->has("JRCS_init_obj_scale"))
    {
        init_obj_scale_ = config_->getFloat("JRCS_init_obj_scale");
    }else init_obj_scale_ = 1.0;
    return true;
}


void JRCSBase::reset_iteration()
{
    iter_count_ = 0;
}

void JRCSBase::input(
        const MatPtrLst& vv,
        const MatPtrLst& vn,
        const CMatPtrLst& vc,
        const LCMatPtrLst& vl
        )
{
    vvs_ptrlst_ = vv;
    vns_ptrlst_ = vn;
    vcs_ptrlst_ = vc;
    vls_ptrlst_ = vl;
}

void JRCSBase::input_with_label(
        const MatPtrLst& vv,
        const MatPtrLst& vn,
        const CMatPtrLst& vc,
        const LCMatPtrLst& vlc,
        const LMatPtrLst& vl
        )
{
    vvs_ptrlst_ = vv;
    vns_ptrlst_ = vn;
    vcs_ptrlst_ = vc;
    vls_ptrlst_ = vlc;
    vll_ptrlst_ = vl;
}

void JRCSBase::resetw(
        const MatPtrLst& wv,
        const MatPtrLst& wn,
        const CMatPtrLst& wc
        )
{
    wvs_ptrlst_ = wv;
    wns_ptrlst_ = wn;
    wcs_ptrlst_ = wc;
    int k = wvs_ptrlst_.front()->n_cols;

    MatPtrLst::iterator wviter = wvs_ptrlst_.begin();
    MatPtrLst::iterator vviter = vvs_ptrlst_.begin();
    MatPtrLst::iterator wniter = wns_ptrlst_.begin();
    MatPtrLst::iterator vniter = vns_ptrlst_.begin();
    CMatPtrLst::iterator wciter = wcs_ptrlst_.begin();
    CMatPtrLst::iterator vciter = vcs_ptrlst_.begin();
    for( int idx = 0; idx < wvs_ptrlst_.size() ; ++ idx )
    {
        int range = ( *vviter )->n_cols - 1;
        arma::uvec index = arma::randi<arma::uvec>( k, arma::distr_param( 0 , range ) );

        **wviter = (**vviter).cols(index);
        **wniter = (**vniter).cols(index);
        **wciter = (**vciter).cols(index);

        ++wviter;
        ++vviter;
        ++wniter;
        ++vniter;
        ++wciter;
        ++vciter;
    }
}

void JRCSBase::reset_objw(const std::vector<float>& pvec)
{
    std::cerr<<"JRCSBase::reset_objw"<<std::endl;
    obj_num_ = pvec.size();
    obj_prob_= arma::fmat(pvec);
    obj_prob_ = arma::sort(obj_prob_);
}

void JRCSBase::initx(
        const MatPtr& xv,
        const MatPtr& xn,
        const CMatPtr& xc
        )
{
    xv_ptr_ = xv;
    xn_ptr_ = xn;
    xc_ptr_ = xc;
    xtv_ = *xv_ptr_;
    xtn_ = *xn_ptr_;
    xtc_ = *xc_ptr_;

    int k = xv_ptr_->n_cols;
    if(verbose_>0)std::cerr<<"k:"<<k<<std::endl;

    init_alpha_ = false;
    if(!vll_ptrlst_.empty())
    {
        if(!init_)init_.reset(new JRCSInitBase());
        init_->configure(config_);
        if(init_->init_with_label(
                    k,
                    vvs_ptrlst_,
                    vns_ptrlst_,
                    vcs_ptrlst_,
                    vls_ptrlst_,
                    vll_ptrlst_,
                    verbose_)
           )init_alpha_=true;
    }
    if(verbose_>0)std::cerr<<"manually obj_prob:"<<std::endl;
    if(verbose_>0)std::cerr<<obj_prob_<<std::endl;
    if(init_alpha_)
    {
        init_->getObjProb(obj_prob_);
        obj_num_ = obj_prob_.size();
    }
    if(verbose_>0)
    {
        std::cerr<<"manually obj_prob is ignored"<<std::endl;
        std::cerr<<"new obj_prob:"<<std::endl;
        std::cerr<<obj_prob_<<std::endl;
    }
    int r_k = k;
    float* pxv = (float*)xtv_.memptr();
    float* pxn = (float*)xtn_.memptr();
    uint8_t* pxc = (uint8_t*)xtc_.memptr();
    int N = obj_prob_.size();
    obj_pos_ = arma::fmat(3,N,arma::fill::zeros);
    arma::frowvec z = arma::linspace<arma::frowvec>(float(-N/2),float(N/2),N);
    obj_pos_.row(2) = z;
    max_obj_radius_ = 0.0;
    if(verbose_>0){
        std::cerr<<"obj_pos_:"<<std::endl;
        std::cerr<<obj_pos_<<std::endl;
    }
    obj_label_ = arma::uvec(xv_ptr_->n_cols,arma::fill::zeros);
    arma::uword* pxl = (arma::uword*)obj_label_.memptr();
    for(int obj_idx = 0 ; obj_idx < obj_prob_.size() ; ++ obj_idx )
    {
        int obj_size = int(float(k)*float(obj_prob_(obj_idx)));
        obj_size = std::max(9,obj_size);
        obj_size = std::min(r_k,obj_size);
        if( obj_idx == ( obj_prob_.size() - 1 ) ) obj_size = r_k;
        objv_ptrlst_.emplace_back(new arma::fmat(pxv,3,obj_size,false,true));
        objn_ptrlst_.emplace_back(new arma::fmat(pxn,3,obj_size,false,true));
        objc_ptrlst_.emplace_back(new arma::Mat<uint8_t>(pxc,3,obj_size,false,true));
        arma::uvec lbl(pxl,obj_size,false,true);
        lbl.fill(arma::uword(obj_idx+1));
        pxv += 3*obj_size;
        pxn += 3*obj_size;
        pxc += 3*obj_size;
        pxl += obj_size;
        r_k -= obj_size;
        arma::fvec pos = obj_pos_.col(obj_idx);
        reset_obj_vn(1.0,pos,(*objv_ptrlst_.back()),(*objn_ptrlst_.back()));
        reset_obj_c((*objc_ptrlst_.back()));
    }
    *xv_ptr_ = xtv_;
    *xn_ptr_ = xtn_;
    *xc_ptr_ = xtc_;
    if(verbose_)std::cerr<<"Done initx_"<<std::endl;
}

void JRCSBase::reset_obj_vn(
        float radius,
        arma::fvec& pos,
        arma::fmat& ov,
        arma::fmat& on
        )
{
    rand_sphere(ov);
    on = ov;
    ov *= radius;
    if( radius > max_obj_radius_ ) max_obj_radius_ = radius;
    ov.each_col() += pos;
}

void JRCSBase::rand_sphere(
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
}

void JRCSBase::reset_obj_c(
        arma::Mat<uint8_t>& oc
        )
{
    oc.fill(128);
}

void JRCSBase::reset_rt()
{
    if(vvs_ptrlst_.empty())
    {
        throw std::logic_error("JRCSBase::reset_rt: no input");
    }
    std::cerr<<"JRCSBase::reset_rt()"<<std::endl;
    rt_lst_.resize(vvs_ptrlst_.size());
    TsLst::iterator iter;
    for(iter=rt_lst_.begin();iter!=rt_lst_.end();++iter)
    {
        Ts& rt = *iter ;
        rt.resize(obj_num_);
        #pragma omp parallel for
        for( int i = 0 ; i < obj_num_ ; ++i )
        {
            arma::fmat R(rt[i].R,3,3,false,true);
            arma::fvec t(rt[i].t,3,false,true);
            R.fill(arma::fill::eye);
            arma::fvec randw(2*obj_num_,arma::fill::randu);
            randw = arma::normalise(randw);
            arma::fvec randt(3,arma::fill::randn);
            arma::fmat pos0 = obj_pos_.each_col() + randt*max_obj_radius_;
            arma::fmat pos1 = obj_pos_.each_col() - randt*max_obj_radius_;
            arma::fmat pos_base = arma::join_rows(pos0,pos1);
            t = pos_base*randw;
//            std::cerr<<"t for obj["<<i<<"]:"<<t.t()<<std::endl;
        }
    }
}

int JRCSBase::evaluate_k(const arma::uvec& sizes)
{
    arma::uword k = arma::median(sizes);
    return ( k/2 + 7 ) ;
}

int JRCSBase::evaluate_k()
{
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
    int k =  evaluate_k(k_lst);//median size but at least five;
    return k;
}

void JRCSBase::computeOnce()
{
    xv_sum_.fill(0.0);
    xn_sum_.fill(0.0);
    xc_sum_.fill(0.0);
    var_sum.fill(0.0);
    alpha_sum.fill(0.0);
    alpha_sumij.fill(0.0);

    //reset transformed latent center
    if(verbose_>0)std::cerr<<"reset transformed latent color"<<std::endl;
    xtc_ = *xc_ptr_;

    if(smooth_enabled_)computeCompatibility(mu_);

    if(verbose_>1)
    {
        std::stringstream muname;
        muname.str("");
        muname<<debug_path_<<"mu_"<<iter_count_<<".fmat.arma";
        mu_.save(muname.str(),arma::raw_ascii);
    }

    for(int idx=0;idx<vvs_ptrlst_.size();++idx)
    {
        //update rate for this frame
        float rate0 = float(idx) / float(idx+1);
        float rate1 = 1.0 - rate0;
        if(verbose_>0)std::cerr<<"reset transformed latent center"<<std::endl;
        xtv_ = *xv_ptr_;
        xtn_ = *xn_ptr_;

        arma::fmat& vv_ = *vvs_ptrlst_[idx];
        arma::fmat& vn_ = *vns_ptrlst_[idx];
        arma::Mat<uint8_t>& vc_ = *vcs_ptrlst_[idx];
        arma::mat& alpha = *alpha_ptrlst_[idx];
        Ts& rt = rt_lst_[idx];

        if(verbose_>0)std::cerr<<"step-E"<<std::endl;
        if(verbose_>0)std::cerr<<"transform object"<<std::endl;
        #pragma omp parallel for
        for(int o = 0 ; o < obj_num_ ; ++o )
        {
            arma::fmat R(rt[o].R,3,3,false,true);
            arma::fvec t(rt[o].t,3,false,true);
            arma::fmat& objv = *objv_ptrlst_[o];
            arma::fmat& objn = *objn_ptrlst_[o];
            objv = R*objv;
            objv.each_col() += t;
            objn = R*objn;
        }

        arma::mat tmpxc = arma::conv_to<arma::mat>::from(xtc_);
        arma::mat tmpvc = arma::conv_to<arma::mat>::from(vc_);

        if(verbose_>0)std::cerr<<"calculate alpha"<<std::endl;
        if( (iter_count_>0) || (!init_alpha_) )
        {
            #pragma omp parallel for
            for(int r = 0 ; r < alpha.n_rows ; ++r )
            {
                arma::mat tmpv = arma::conv_to<arma::mat>::from(xtv_.each_col() - vv_.col(r)); //after the objv is transformed the xtv is transformed
                alpha.row(r)  = arma::sum(arma::square(tmpv));
            }

            alpha.each_row() %= (-0.5*x_invvar_);
            alpha = arma::trunc_exp(alpha);
            alpha.each_row() %= arma::pow(x_invvar_,1.5);
            alpha.each_row() %= x_p_;
        }else{
            if(verbose_>0)std::cerr<<"using init alpha"<<std::endl;
        }

        //normalise alpha
        arma::vec alpha_rowsum = arma::sum(alpha,1) + beta_;
        alpha.each_col() /= alpha_rowsum;
        alpha_rowsum = arma::sum(alpha,1);

        if(verbose_>1)
        {
            std::stringstream alphaname;
            alphaname.str("");
            alphaname<<debug_path_<<"alpha_"<<idx<<"_iter_"<<iter_count_<<".fmat.arma";
            alpha.save(alphaname.str(),arma::raw_ascii);
        }
        //smoothing alpha
        if(smooth_enabled_ && iter_count_ > max_init_iter_)
        {
            if(verbose_>0)std::cerr<<"smoothing alpha"<<std::endl;

            DenseCRF3D crf(vv_,vn_,arma::conv_to<arma::fmat>::from(tmpvc),xv_ptr_->n_cols);
            arma::mat unary = arma::conv_to<arma::mat>::from(-1.0*arma::log(alpha));

            crf.setUnaryEnergy(unary.t());
            arma::fvec sxyz = { 0.05 , 0.05 , 0.05 } ;
            crf.addPairwiseGaussian(sxyz,new MatrixCompatibility(mu_));
            arma::fvec srgb = { 10 , 10 , 10 };
            arma::fvec snxyz = { 0.05 , 0.05, 0.05 };
            crf.addPairwiseBilateral(sxyz,snxyz,srgb,new MatrixCompatibility(mu_));

            if(verbose_>0)std::cerr<<"start smoothing"<<std::endl;
            arma::mat Q = crf.startInference();
            arma::mat t1,t2;
            if(verbose_>0)std::cerr<<"kl = "<<crf.klDivergence(Q)<<std::endl;
            for( int it=0; it<max_smooth_iter_; it++ ) {
                crf.stepInference( Q, t1, t2 );
                if(verbose_>0)std::cerr<<"kl = "<<crf.klDivergence(Q)<<std::endl;
            }
            alpha = arma::conv_to<arma::mat>::from(Q.t());
            alpha_rowsum = ( 1.0 + beta_ ) * arma::sum(alpha,1);
            alpha.each_col() /= alpha_rowsum;
            alpha_rowsum = arma::sum(alpha,1);
            if(verbose_>1)
            {
                std::stringstream alphaname;
                alphaname.str("");
                alphaname<<debug_path_<<"alpha_"<<idx<<"_iter_"<<iter_count_<<"_smooth.fmat.arma";
                alpha.save(alphaname.str(),arma::raw_ascii);
            }
        }
        //update RT
        //#1 calculate weighted point cloud
        if(verbose_>0)std::cerr<<"calculating the weighted point cloud"<<std::endl;
        arma::fmat& wv = *wvs_ptrlst_[idx];
        arma::fmat& wn = *wns_ptrlst_[idx];
        arma::fmat wc = arma::conv_to<arma::fmat>::from(*wcs_ptrlst_[idx]);
        arma::rowvec alpha_colsum = arma::sum( alpha );
        arma::rowvec alpha_median = arma::median( alpha );
        arma::mat trunc_alpha = alpha;

        #pragma omp parallel for
        for(int c=0;c<alpha.n_cols;++c)
        {
            arma::vec col = trunc_alpha.col(c);
            col( col < alpha_median(c) ).fill(0.0);
            trunc_alpha.col(c) = col;
        }

        trunc_alpha += std::numeric_limits<double>::epsilon(); //add eps for numerical stability

        arma::rowvec trunc_alpha_colsum = arma::sum(trunc_alpha);

        arma::frowvec square_lambda = arma::conv_to<arma::frowvec>::from(x_invvar_ % trunc_alpha_colsum);
        arma::frowvec p(square_lambda.n_cols,arma::fill::ones);
        arma::frowvec square_norm_lambda = square_lambda / arma::accu(square_lambda);
        p -= square_norm_lambda;

        wv = vv_*arma::conv_to<arma::fmat>::from(trunc_alpha);

        wn = vn_*arma::conv_to<arma::fmat>::from(trunc_alpha);
        wc = arma::conv_to<arma::fmat>::from(vc_)*arma::conv_to<arma::fmat>::from(trunc_alpha);

        #pragma omp parallel for
        for(int c=0;c<alpha.n_cols;++c)
        {
            if( 0 != trunc_alpha_colsum(c) )
            {
                wv.col(c) /= trunc_alpha_colsum(c);
                wn.col(c) /= trunc_alpha_colsum(c);
                wc.col(c) /= trunc_alpha_colsum(c);
            }
        }

        wn = arma::normalise( wn );
        *wcs_ptrlst_[idx] = arma::conv_to<arma::Mat<uint8_t>>::from(wc);

//        arma::fmat closest_v = vv_.cols(closest_i);
//        arma::fmat closest_n = vn_.cols(closest_i);
//        arma::Mat<uint8_t> closest_c8 = vc_.cols(closest_i);
//        arma::fmat closest_c = arma::conv_to<arma::fmat>::from(closest_c8);

        if(verbose_>0)std::cerr<<"calculating R & t"<<std::endl;
        if(verbose_>0)std::cerr<<"rate0:"<<rate0<<","<<"rate1:"<<rate1<<std::endl;
        #pragma omp parallel for
        for(int o = 0 ; o < obj_num_ ; ++o )
        {
            arma::fmat A;
            arma::fmat U,V;
            arma::fvec s;
            arma::fmat R(rt[o].R,3,3,false,true);
            arma::fvec t(rt[o].t,3,false,true);
            arma::fmat dR;
            arma::fvec dt;
            arma::fmat objv = *objv_ptrlst_[o];
            arma::uvec oidx = arma::find(obj_label_==(o+1));
            if(verbose_>1)std::cerr<<"(s,e)=("<<oidx.min()<<","<<oidx.max()<<")"<<std::endl;
            arma::fmat v;
            v = wv.cols(oidx);
            arma::fmat cv = v.each_col() - arma::mean(v,1);
            objv.each_col() -= arma::mean(objv,1);
            objv.each_row() %= p.cols(oidx) % square_lambda.cols(oidx) ;
            A = cv*objv.t();
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
                    arma::fmat ddt = v - dR*(*objv_ptrlst_[o]);
                    ddt.each_row() %= square_norm_lambda.cols(oidx);
                    dt = arma::sum(ddt,1);
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
                    arma::fmat ddt = v - dR*(*objv_ptrlst_[o]);
                    ddt.each_row() %= square_norm_lambda.cols(oidx);
                    dt = arma::sum(ddt,1);
                }
            }
            }

            //updating objv
            //not needed since the transformed obj will be reset
            /*
            *objv_ptrlst_[o] = dR*(*objv_ptrlst_[o]);
            (*objv_ptrlst_[o]).each_col() += dt;
            *objn_ptrlst_[o] = dR*(*objn_ptrlst_[o]);
            */

            //updating R T
            R = dR*R;
            t = dR*t + dt;

            //accumulate for updating X
            arma::fmat tv = vv_.each_col() - t;
            tv = R.i() * tv;
            arma::fmat twv = tv*arma::conv_to<arma::fmat>::from(trunc_alpha);
            xv_sum_.cols(oidx) +=  twv.cols(oidx);
            xn_sum_.cols(oidx) = rate0*xn_sum_.cols(oidx)+rate1*R.i()*wn.cols(oidx);
            xc_sum_.cols(oidx) = rate0*xc_sum_.cols(oidx)+rate1*wc.cols(oidx);
        }
        //update var
        alpha_sum += trunc_alpha_colsum;
        arma::fmat alpha_2(alpha.n_rows,alpha.n_cols);
        #pragma omp parallel for
        for(int r=0;r<alpha_2.n_rows;++r)
        {
            alpha_2.row(r) = arma::sum(arma::square(xtv_.each_col() - vv_.col(r)));
        }
        arma::rowvec tmpvar = arma::sum(alpha_2%alpha);
        var_sum += tmpvar;
        alpha_sumij += alpha_colsum;
        QCoreApplication::processEvents();
    }

    if(verbose_>0)std::cerr<<"Updating X:"<<std::endl;
//    float N =  vvs_ptrlst_.size();
    assert(xv_sum_.is_finite());
    xv_sum_.each_row() /= arma::conv_to<arma::frowvec>::from(alpha_sum);
    *xv_ptr_ = xv_sum_;
//    assert((*xv_ptr_).has_inf()||(*xv_ptr_).has_nan());
    //fix the x center position
//    std::cerr<<"xv:"<<xv_ptr_->n_cols<<std::endl;
    #pragma omp parallel for
    for(int o = 0 ; o < obj_num_ ; ++o )
    {
        arma::uvec oidx = arma::find(obj_label_==(o+1));
        arma::fmat newxv = xv_ptr_->cols(oidx);
        arma::fvec adjustt =  obj_pos_.col(o) - arma::mean(newxv,1);
        xv_ptr_->cols(oidx) = newxv.each_col() + adjustt;
        for(int idx=0;idx<vvs_ptrlst_.size();++idx)
        {
            arma::fvec t(rt_lst_[idx][o].t,3,false,true);
            t += adjustt;
        }
    }

    *xn_ptr_ = xn_sum_;
    *xn_ptr_ = arma::normalise(*xn_ptr_);
    *xc_ptr_ = arma::conv_to<arma::Mat<uint8_t>>::from( xc_sum_ );

    x_invvar_ = ( (3.0*alpha_sum ) / ( var_sum + 1e-6 ) );//restore reciprocal fractions of variation
    float mu = arma::accu(alpha_sum);
    mu *= ( 1.0 + beta_ );
    x_p_ = alpha_sum;
    if( mu != 0)x_p_ /= mu;
}

void JRCSBase::obj_only(arma::mat& mu)
{
    #pragma omp parallel for
    for(int o = 0 ; o < obj_num_ ; ++o )
    {
        arma::uvec oidx = arma::find(obj_label_==(o+1));
        arma::mat m;
        m = mu(oidx,oidx);
        m.fill( -smooth_w_ );
        mu(oidx,oidx) = m ;
    }
}

void JRCSBase::obj_point_dist(arma::mat& mu)
{
    #pragma omp parallel for
    for(int o = 0 ; o < obj_num_ ; ++o )
    {
        arma::uvec oidx = arma::find(obj_label_==(o+1));
        arma::fmat xv = xv_ptr_->cols(oidx);
        arma::mat m;
        m = mu(oidx,oidx);
        for(int i=0 ; i < xv.n_cols ; ++i )
        {
            for(int j = i ; j < xv.n_cols ; ++j )
            {
                arma::fvec dx = xv.col(i) - xv.col(j);
                m(i,j) = -smooth_w_ * std::exp( - arma::sum( arma::square( dx / 0.5 ) ));
            }
        }
        mu(oidx,oidx) = m ;
    }
}

void JRCSBase::computeCompatibility(arma::mat& mu)
{
    if(verbose_>0)std::cerr<<"computeCompatibility"<<std::endl;
    if(verbose_>0)std::cerr<<"smooth_weight:"<<smooth_w_<<std::endl;
    mu = arma::mat(xv_ptr_->n_cols,xv_ptr_->n_cols,arma::fill::zeros);
    switch(mu_type_)
    {
    case ObjOnly:
        obj_only(mu);
        break;
    case ObjPointDist:
        obj_point_dist(mu);
        break;
    }
}

void JRCSBase::reset_alpha()
{
    if(init_alpha_)
    {
       assert(init_&&(init_.use_count()>0));
       init_->getAlpha(alpha_ptrlst_);
    }else
    {
        if(verbose_>0)std::cerr<<"allocating alpha"<<std::endl;
        arma::fmat& xv_ = *xv_ptr_;
        int idx=0;
        while( idx < vvs_ptrlst_.size() )
        {
            if(idx>=alpha_ptrlst_.size())alpha_ptrlst_.emplace_back(new arma::mat(vvs_ptrlst_[idx]->n_cols,xv_.n_cols));
            else if((alpha_ptrlst_[idx]->n_rows!=vvs_ptrlst_[idx]->n_cols)||(alpha_ptrlst_[idx]->n_cols!=xv_.n_cols))
            alpha_ptrlst_[idx].reset(new arma::mat(vvs_ptrlst_[idx]->n_cols,xv_.n_cols));
            ++idx;
        }
        if(verbose_>0)std::cerr<<"done allocating alpha"<<std::endl;
    }
}

void JRCSBase::reset_prob()
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

    float maxvar = arma::accu(arma::square(maxAllXYZ-minAllXYZ));
    (*xv_ptr_) *= (0.5*std::sqrt(maxvar));
    x_invvar_ = arma::rowvec(xv_ptr_->n_cols);
    x_invvar_.fill(1.0/maxvar);

    xv_sum_ = arma::fmat(xv_ptr_->n_rows,xv_ptr_->n_cols,arma::fill::zeros);
    xn_sum_ = arma::fmat(xn_ptr_->n_rows,xn_ptr_->n_cols,arma::fill::zeros);
    xc_sum_ = arma::fmat(xc_ptr_->n_rows,xc_ptr_->n_cols,arma::fill::zeros);

    x_p_ = arma::rowvec(xv_ptr_->n_cols);
    x_p_.fill(1.0/float(xv_ptr_->n_cols));

    var_sum = arma::rowvec(xv_ptr_->n_cols,arma::fill::zeros);
    alpha_sum = arma::rowvec(xv_ptr_->n_cols,arma::fill::zeros);
    alpha_sumij = arma::rowvec(xv_ptr_->n_cols,arma::fill::zeros);

    float h;
    h = 2.0 / arma::mean(x_invvar_);
    beta_ = 0.1/(h*(0.1+1));
    if(verbose_>0)std::cerr<<"done probability"<<std::endl;
}

void JRCSBase::update_color_label()
{
    if(verbose_>0)std::cerr<<"updating color label"<<std::endl;
    for(int idx=0;idx<vvs_ptrlst_.size();++idx)
    {
        arma::mat& alpha = *alpha_ptrlst_[idx];
        arma::mat obj_p(alpha.n_rows,obj_num_);
        arma::Col<uint32_t>& vl = *vls_ptrlst_[idx];
        #pragma omp parallel for
        for(int o = 0 ; o < obj_num_ ; ++o )
        {
            arma::uvec oidx = arma::find(obj_label_==(o+1));
            arma::mat sub_alpha = alpha.cols(oidx);
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

void JRCSBase::get_label(std::vector<arma::uvec>&lbl)
{
    if(verbose_)std::cerr<<"JRCSBase::get_label"<<std::endl;
    if(lbl.size()!=vvs_ptrlst_.size()){
        lbl.resize(vvs_ptrlst_.size());
    }
    for(int idx=0;idx<vvs_ptrlst_.size();++idx)
    {
        arma::mat& alpha = *alpha_ptrlst_[idx];
        arma::mat obj_p(alpha.n_rows,obj_num_);
        arma::Col<uint32_t>& vl = *vls_ptrlst_[idx];
        #pragma omp parallel for
        for(int o = 0 ; o < obj_num_ ; ++o )
        {
            arma::uvec oidx = arma::find(obj_label_==(o+1));
            arma::mat sub_alpha = alpha.cols(oidx);
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
        lbl[idx] = label;
    }
}

typedef struct{
    int idx;
    arma::fvec coord;
}Xorder;

arma::fvec order_eps;

bool ordering(const Xorder& x1,const Xorder& x2)
{
    if( std::fabs( x1.coord(2) - x2.coord(2) ) < order_eps(2) ) // consider equal
    {
        if( std::fabs( x1.coord(1) - x2.coord(1) ) < order_eps(1) )
        {
            return x1.coord(0) < x2.coord(0);
        }else return x1.coord(1) < x2.coord(1);
    }
    else return x1.coord(2) < x2.coord(2);
}

void JRCSBase::get_order(std::vector<arma::uvec>& orders)
{
    if(verbose_)std::cerr<<"JRCSBase::get_order"<<std::endl;
    if(orders.size()!=(vvs_ptrlst_.size()+1)){
        orders.resize(vvs_ptrlst_.size()+1);
    }
    arma::fvec maxX = arma::max(*xv_ptr_,1);
    arma::fvec minX = arma::min(*xv_ptr_,1);
    order_eps = ( maxX - minX ) / 20.0;
    std::vector<Xorder> Xscore(xv_ptr_->n_cols);
    #pragma omp parallel for
    for(int i = 0 ; i < Xscore.size() ; ++i )
    {
        Xscore[i].coord = xv_ptr_->col(i);
        Xscore[i].idx = i;
    }
    std::sort(Xscore.begin(),Xscore.end(),ordering);
    arma::uvec Xorder(Xscore.size());
    #pragma omp parallel for
    for(arma::uword i = 0 ; i < Xscore.size() ; ++i )
    {
        Xorder(Xscore[i].idx) = i;
    }
    orders[0] = Xorder;
    for(int idx=0;idx<vvs_ptrlst_.size();++idx)
    {
        arma::mat& alpha = *alpha_ptrlst_[idx];
        arma::uvec order(alpha.n_rows);
        #pragma omp parallel for
        for(int r = 0 ; r < alpha.n_rows ; ++r )
        {
            arma::uword l;
            arma::rowvec point_prob = alpha.row(r);
            point_prob.max(l);
            order(r) = Xorder(l);
        }
        orders[idx+1] = order;
    }
}

bool JRCSBase::isEnd()
{
    if(iter_count_>=max_iter_)return true;
    return false;
}

}
