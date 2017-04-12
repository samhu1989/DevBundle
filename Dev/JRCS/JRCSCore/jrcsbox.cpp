#include "jrcsbox.h"
#include <QThread>
namespace JRCS{

bool JRCSBox::update_cube_;
std::vector<JRCSBox::Cube::PtrLst> JRCSBox::cube_ptrlsts_;

void JRCSBox::set_boxes(std::vector<Cube::PtrLst>& cube_ptrlsts)
{
    cube_ptrlsts_ = cube_ptrlsts;
}

void JRCSBox::add_boxes(std::vector<Cube::PtrLst>& cube_ptrlsts)
{
    update_cube_ = true;
    cube_ptrlsts_ = cube_ptrlsts;
}

JRCSBox::JRCSBox()
{
    ;
}

void JRCSBox::initx(
        const MatPtr& xv,
        const MatPtr& xn,
        const CMatPtr& xc
        )
{
    if(verbose_)std::cerr<<"Start JRCSBox::initx"<<std::endl;
    xv_ptr_ = xv;
    xn_ptr_ = xn;
    xc_ptr_ = xc;
    xtv_ = *xv_ptr_;
    xtn_ = *xn_ptr_;
    xtc_ = *xc_ptr_;

    int k = xv_ptr_->n_cols;
    if(verbose_>0)std::cerr<<"k:"<<k<<std::endl;
    if(verbose_>0)std::cerr<<"obj_num:"<<obj_num_<<std::endl;
    init_from_boxes();

    if(verbose_>0)std::cerr<<"obj_prob:"<<obj_prob_.t()<<std::endl;

    int r_k = k;
    float* pxv = (float*)xtv_.memptr();
    float* pxn = (float*)xtn_.memptr();
    uint8_t* pxc = (uint8_t*)xtc_.memptr();

    int N = obj_num_;
    obj_pos_ = arma::fmat(3,N,arma::fill::zeros);
    arma::frowvec z = arma::linspace<arma::frowvec>(float(-N/2),float(N/2),N);
    obj_pos_.row(2) = z;
    max_obj_radius_ = 0.0;
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
        if( obj_idx == ( obj_prob_.size() - 1 ) ) obj_size = r_k;
        if( pxr > pxr_s )*pxr = (*(pxr - 1))+1;
        *(pxr+1) = *pxr + obj_size - 1;
        objv_ptrlst_.emplace_back(new arma::fmat(pxv,3,obj_size,false,true));
        objn_ptrlst_.emplace_back(new arma::fmat(pxn,3,obj_size,false,true));
        objc_ptrlst_.emplace_back(new arma::Mat<uint8_t>(pxc,3,obj_size,false,true));
        pxv += 3*obj_size;
        pxn += 3*obj_size;
        pxc += 3*obj_size;
        r_k -= obj_size;
        pxr += 2;
        arma::fvec pos = obj_pos_.col(obj_idx);
        reset_obj_vn(0.5,pos,(*objv_ptrlst_.back()),(*objn_ptrlst_.back()));
        reset_obj_c((*objc_ptrlst_.back()));
    }
    *xv_ptr_ = xtv_;
    *xn_ptr_ = xtn_;
    *xc_ptr_ = xtc_;
    if(verbose_)std::cerr<<"Done initx_"<<std::endl;
}

void JRCSBox::update_from_cube(void)
{
    std::vector<Cube::PtrLst>::iterator iter;
    MatPtrLst::iterator vviter = vvs_ptrlst_.begin();
    DMatPtrLst::iterator piter = inbox_prob_lsts_.begin();
    for(iter=cube_ptrlsts_.begin();iter != cube_ptrlsts_.end();++iter)
    {
        if( iter->size() > 0 )
        {
            if(!(*piter))
            {
                if(verbose_)std::cerr<<"updating inbox prob"<<std::endl;
                init_obj_prob(*iter,*vviter,*piter);
            }
        }
        ++vviter;
        if(vviter==vvs_ptrlst_.end())break;
        ++piter;
        if(piter==inbox_prob_lsts_.end())break;
    }
    update_cube_ = false;
}

void JRCSBox::reset_objw(const std::vector<float>&)
{
    if(verbose_)std::cerr<<"JRCSBox::reset_objw"<<std::endl;
    obj_num_ = 0;
    std::vector<Cube::PtrLst>::iterator iter;
    for(iter=cube_ptrlsts_.begin();iter!=cube_ptrlsts_.end();++iter)
    {
        int max_label = 0;
        for(Cube::PtrLst::iterator cube_iter=iter->begin();cube_iter!=iter->end();++cube_iter)
        {
            if(max_label<(*cube_iter)->label_)
            {
                max_label = (*cube_iter)->label_;
            }
        }
        if( obj_num_ < max_label )
        {
            obj_num_ = max_label;
        }
    }
    if(verbose_)std::cerr<<"JRCSBox::reset_objw:obj_num_="<<obj_num_<<std::endl;
}

void JRCSBox::update_color_label()
{
    if(verbose_>0)std::cerr<<"JRCSBox::updating color label"<<std::endl;
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
        Cube::colorByLabel((uint32_t*)vl.memptr(),vl.size(),label);
    }
}

void JRCSBox::reset_alpha()
{
    if(verbose_)std::cerr<<"JRCSBox::reset_alpha():allocating alpha"<<std::endl;
    arma::fmat& xv_ = *xv_ptr_;
    int idx=0;
    while( idx < vvs_ptrlst_.size() )
    {
        if(idx>=alpha_ptrlst_.size())alpha_ptrlst_.emplace_back(new arma::mat(vvs_ptrlst_[idx]->n_cols,xv_.n_cols));
        else if((alpha_ptrlst_[idx]->n_rows!=vvs_ptrlst_[idx]->n_cols)||(alpha_ptrlst_[idx]->n_cols!=xv_.n_cols))
        alpha_ptrlst_[idx].reset(new arma::mat(vvs_ptrlst_[idx]->n_cols,xv_.n_cols));
        ++idx;
    }
    if(verbose_)std::cerr<<"JRCSBox::reset_alpha():done allocating alpha"<<std::endl;
}

void JRCSBox::init_from_boxes()
{
    if(verbose_)std::cerr<<"JRCSBox::init_from_boxes()"<<std::endl;
    std::vector<Cube::PtrLst>::iterator iter;
    MatPtrLst::iterator vviter = vvs_ptrlst_.begin();
    CMatPtrLst::iterator vciter = vcs_ptrlst_.begin();
    inbox_prob_lsts_.resize(vvs_ptrlst_.size());
    DMatPtrLst::iterator piter = inbox_prob_lsts_.begin();
    obj_prob_.clear();
    if(verbose_)std::cerr<<"init inbox prob"<<std::endl;
    for(iter=cube_ptrlsts_.begin();iter!=cube_ptrlsts_.end();++iter)
    {
        if( iter->size() > 0 )
        {
            if(obj_prob_.empty())
            {
                obj_prob_ = obj_prob_from_boxes(*iter,*vviter);
                if(obj_prob_.size()==obj_num_)
                {
                    init_color_gmm(*iter,*vviter,*vciter,color_gmm_lsts_);
                }else obj_prob_.clear();
            }
            init_obj_prob(*iter,*vviter,*piter);
        }
        ++vviter;
        if(vviter==vvs_ptrlst_.end())break;
        ++vciter;
        if(vciter==vcs_ptrlst_.end())break;
        ++piter;
        if(piter==inbox_prob_lsts_.end())break;
    }
    if(verbose_)std::cerr<<"init color prob"<<std::endl;
    color_prob_lsts_.resize(vvs_ptrlst_.size());
    vciter = vcs_ptrlst_.begin();
    for(piter = color_prob_lsts_.begin();piter!=color_prob_lsts_.end();++piter)
    {
        init_color_prob(*vciter,*piter);
        ++vciter;
        if(vciter==vcs_ptrlst_.end())break;
    }
    if(verbose_)std::cerr<<"JRCSBox::init_from_boxes()[END]"<<std::endl;
}

arma::fvec JRCSBox::obj_prob_from_boxes(const Cube::PtrLst& cubes,const MatPtr& vv)
{
    std::vector<float> prob_vec;
    prob_vec.reserve(cubes.size());
    uint32_t idx = 0;
    for( Cube::PtrLst::const_iterator iter = cubes.cbegin() ; iter!=cubes.cend() ; ++iter )
    {
        Cube& cube = **iter;
        if( prob_vec.size() < cube.label_ )prob_vec.push_back(arma::accu( arma::trunc_exp(-cube.get_dist2_box(*vv)) ));
        else prob_vec[cube.label_-1] += arma::accu( arma::trunc_exp(-cube.get_dist2_box(*vv)) );
        ++idx;
    }
    arma::fvec prob(prob_vec);
    return prob / arma::accu(prob);
}

void JRCSBox::init_color_prob(const CMatPtr& c,DMatPtr& prob)
{
    prob.reset(new arma::mat(c->n_cols,color_gmm_lsts_.size()));
    arma::mat data = arma::conv_to<arma::mat>::from(*c);
    #pragma omp parallel for
    for(int i=0;i<prob->n_cols;++i)
    {
        arma::rowvec p = arma::trunc_exp(color_gmm_lsts_[i]->log_p(data));
        prob->col(i) = p.t();
    }
    prob->each_col() /= arma::sum(*prob,1);
}

void JRCSBox::init_color_gmm(
        const Cube::PtrLst& cubes,
        const MatPtr& vv,
        const CMatPtr& vc,
        GMMPtrLst& gmms
        )
{
    gmms.resize(obj_num_);
    for( Cube::PtrLst::const_iterator iter = cubes.cbegin() ; iter!=cubes.cend() ; ++iter )
    {
        Cube& cube = **iter;
        arma::uvec inside = cube.inside(*vv);
        arma::mat data = arma::conv_to<arma::mat>::from( vc->cols(inside) );
        if(!gmms[cube.label_-1])
        {
            gmms[cube.label_-1].reset(new arma::gmm_diag());
            if(verbose_)gmms[cube.label_-1]->learn(data,5,arma::maha_dist,arma::random_spread,0,30,1e-10,true);
            else gmms[cube.label_-1]->learn(data,5,arma::maha_dist,arma::random_spread,0,30,1e-10,false);
        }else{
            if(verbose_)gmms[cube.label_-1]->learn(data,5,arma::maha_dist,arma::keep_existing,0,30,1e-10,true);
            else gmms[cube.label_-1]->learn(data,5,arma::maha_dist,arma::keep_existing,0,30,1e-10,false);
        }
    }
}

void JRCSBox::init_obj_prob(
        const Cube::PtrLst& cubes,
        const MatPtr& vv,
        DMatPtr& prob
        )
{
    prob.reset(new arma::mat(vv->n_cols,obj_num_,arma::fill::zeros));
    int idx = 0;
    for( Cube::PtrLst::const_iterator iter = cubes.cbegin() ; iter!=cubes.cend() ; ++iter )
    {
        Cube& cube = **iter;
        arma::vec p = arma::trunc_exp( - cube.get_dist2_box(*vv) );
        #pragma omp parallel for
        for(int i=0;i<p.size();++i)
        {
            (*prob)(i,cube.label_-1) = std::max(p(i),(*prob)(i,cube.label_-1));
        }
        ++idx;
    }
}

void JRCSBox::prepare_compute(void)
{
    if(verbose_)std::cerr<<"JRCSBox::prepare_compute"<<std::endl;
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
    if(verbose_){
        std::cerr<<"JRCSBox::prepare_compute[END]"<<std::endl;
        debug_inbox_prob();
        QThread::currentThread()->sleep(5);
        debug_color_prob();
        QThread::currentThread()->sleep(5);
        std::cerr<<"JRCSBox::beta_:"<<beta_<<std::endl;
    }
}

void JRCSBox::step_a(int i)
{
    //latent model
    arma::fmat& xv_ = *xv_ptr_;
    arma::fmat& xn_ = *xn_ptr_;
//    arma::fmat& xf_ = *xf_ptr_;

    //input
    arma::fmat& vv_ = *vvs_ptrlst_[i];
    arma::fmat& vn_ = *vns_ptrlst_[i];
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
        arma::mat tmpv = arma::conv_to<arma::mat>::from( xtv_.each_col() - vv_.col(r) );
        arma::rowvec alpha_v = arma::sum(arma::square(tmpv))%(-0.5*xv_invvar_);
        alpha_v = arma::trunc_exp(alpha_v);
        alpha_v %= arma::pow(xv_invvar_,1.5);
        alpha.row(r) = alpha_v;
    }

     //applying terms of init alpha
    if(iter_count_ < max_init_iter_)
    {
        arma::mat& c_alpha  = *color_prob_lsts_[i];
        int o = 0;
        for(int c = 0 ; c < alpha.n_cols ; ++c )
        {
            alpha.col(c) %= c_alpha.col(o);
            if( c == obj_range_[2*o+1] )++o;//reach end of this object
        }
    }

    //applying terms of points in box constraint to alpha
    if(inbox_prob_lsts_[i])
    {
        arma::mat& b_alpha  = *inbox_prob_lsts_[i];
        int o = 0;
        for(int c = 0 ; c < alpha.n_cols ; ++c )
        {
            alpha.col(c) %= b_alpha.col(o);
            if( c == obj_range_[2*o+1] )++o;//reach end of this object
        }
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
    arma::frowvec alpha_colsum = arma::conv_to<arma::frowvec>::from(arma::sum( alpha ));
    calc_weighted(vv_,vn_,vc,alpha,wv,wn,wc);

    if(verbose_>1)std::cerr<<"#2 calculating RT for each object"<<std::endl;
    for(int o = 0 ; o < obj_num_ ; ++o )
    {
        arma::fmat _v = wv.cols(obj_range_[2*o],obj_range_[2*o+1]);
        arma::fmat _n = wn.cols(obj_range_[2*o],obj_range_[2*o+1]);
        arma::fmat R(rt[o].R,3,3,false,true);
        arma::fvec t(rt[o].t,3,false,true);
        arma::uword offset = 3*obj_range_[2*o];
        arma::uword size = obj_range_[2*o+1] - obj_range_[2*o] + 1;
        arma::fmat objv(((float*)xtv_.memptr())+offset,3,size,false,true);
        arma::fmat objn(((float*)xtn_.memptr())+offset,3,size,false,true);
        updateRTforObj(obj_range_[2*o],obj_range_[2*o+1],_v,_n,objv,objn,alpha_colsum,R,t);
    }
    if(verbose_>1)std::cerr<<"#3 done RT for each object"<<std::endl;

    arma::fmat alpha_v2(alpha.n_rows,alpha.n_cols);
    for(int r=0;r<alpha_v2.n_rows;++r)
    {
        alpha_v2.row(r) = arma::sum(arma::square(xtv_.each_col() - vv_.col(r)));
    }
    vvar_.row(i) = arma::sum(alpha_v2%alpha);
}

void JRCSBox::step_b(void)
{
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

    if(verbose_>1)std::cerr<<"Updating weight of centroids"<<std::endl;
    float mu = arma::accu(alpha_sum);
    mu *= ( 1.0 + beta_ );
    x_p_ = alpha_sum;
    if( mu != 0)x_p_ /= mu;
}

void JRCSBox::finish_steps()
{
    if(update_cube_)update_from_cube();
    JRCSBilateral::finish_steps();
}

void JRCSBox::calc_weighted(
        const arma::fmat&vv,
        const arma::fmat&vn,
        arma::Mat<uint8_t>&vc,
        const arma::mat& alpha,
        arma::fmat&wv,
        arma::fmat&wn,
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
    assert(wv.is_finite());
    //normal is weighted differently
    for(int c=0;c<wn.n_cols;++c)
    {
        arma::mat wvn = arma::conv_to<arma::mat>::from(vn);
        if(!wvn.is_finite())
        {
            arma::uvec idx = find_nonfinite( wvn );
            std::cerr<<idx.size()<<" infinite number of wvn"<<std::endl;
        }
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
    assert(wn.is_finite());

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
            fwc.col(c) /= trunc_alpha_colsum(c);
        }
    }

    if(!wv.is_finite())
    {
        std::cerr<<iter_count_<<":b:!wv.is_finite()"<<std::endl;
    }

    wc = arma::conv_to<arma::Mat<uint8_t>>::from(fwc);
    if(verbose_>1)std::cerr<<"JRCSBilateral::calc_weighted finished"<<std::endl;
}


void JRCSBox::updateRTforObj(
        const arma::uword start,
        const arma::uword end,
        arma::fmat& vv,
        arma::fmat& vn,
        arma::fmat& objv,
        arma::fmat& objn,
        arma::frowvec& alpha_colsum,
        arma::fmat& R,
        arma::fvec& t
        )
{
    arma::fmat A;
    arma::fmat U,V;
    arma::fvec s;
    arma::fmat dR;
    arma::fvec dt;

    arma::fmat cv = vv.each_col() - arma::mean(vv,1);
    if(!vv.is_finite())
    {
        std::cerr<<"JRCSBox::updateRTforObj:!v.is_finite()"<<std::endl;
    }

    arma::frowvec square_lambdav = arma::conv_to<arma::frowvec>::from(xv_invvar_.cols(start,end));
    square_lambdav %= alpha_colsum.cols(start,end);
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
            arma::fmat tmp = vv - dR*objv;
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
            arma::fmat tmp = vv - dR*objv;
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
    objv = dR*objv;
    objv.each_col() += dt;

    objn = dR*objn;

    //updating R T
    R = dR*R;
    t = dR*t + dt;
}

void JRCSBox::debug_inbox_prob()
{
    if(verbose_>0)std::cerr<<"JRCSBox::debug_inbox_prob()"<<std::endl;
    for(int idx=0;idx<vvs_ptrlst_.size();++idx)
    {
        arma::mat& alpha = *alpha_ptrlst_[idx];
        arma::Col<uint32_t>& vl = *vls_ptrlst_[idx];
        if(inbox_prob_lsts_[idx])
        {
            arma::mat& obj_p = *inbox_prob_lsts_[idx];
            arma::uvec label(alpha.n_rows);
            #pragma omp parallel for
            for(int r = 0 ; r < obj_p.n_rows ; ++r )
            {
                arma::uword l;
                arma::rowvec point_prob = obj_p.row(r);
                point_prob.max(l);
                label(r) = l+1;
            }
            Cube::colorByLabel((uint32_t*)vl.memptr(),vl.size(),label);
        }
    }
}

void JRCSBox::debug_color_prob()
{
    if(verbose_>0)std::cerr<<"JRCSBox::debug_color_prob()"<<std::endl;

    for(int idx=0;idx<vvs_ptrlst_.size();++idx)
    {
        arma::mat& alpha = *alpha_ptrlst_[idx];
        arma::Col<uint32_t>& vl = *vls_ptrlst_[idx];
        arma::mat& obj_p  = *color_prob_lsts_[idx];
        arma::uvec label(alpha.n_rows);
        #pragma omp parallel for
        for(int r = 0 ; r < obj_p.n_rows ; ++r )
        {
            arma::uword l;
            arma::rowvec point_prob = obj_p.row(r);
            point_prob.max(l);
            label(r) = l+1;
        }
        Cube::colorByLabel((uint32_t*)vl.memptr(),vl.size(),label);
    }
}

}
