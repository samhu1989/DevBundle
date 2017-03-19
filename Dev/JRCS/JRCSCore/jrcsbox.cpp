#include "jrcsbox.h"
namespace JRCS{

std::vector<JRCSBox::Cube::PtrLst> JRCSBox::cube_ptrlsts_;

void JRCSBox::set_boxes(std::vector<Cube::PtrLst>& cube_ptrlsts)
{
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

    init_from_boxes();

    if(verbose_>0)std::cerr<<obj_prob_<<std::endl;

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
        reset_obj_vn(0.5,pos,(*objv_ptrlst_.back()),(*objn_ptrlst_.back()));
        reset_obj_c((*objc_ptrlst_.back()));
    }
    *xv_ptr_ = xtv_;
    *xn_ptr_ = xtn_;
    *xc_ptr_ = xtc_;
    if(verbose_)std::cerr<<"Done initx_"<<std::endl;
}

void JRCSBox::reset_objw(const std::vector<float>&)
{
    std::cerr<<"JRCSBox::reset_objw"<<std::endl;
    obj_num_ = 0;
    std::vector<Cube::PtrLst>::iterator iter;
    for(iter=cube_ptrlsts_.begin();iter!=cube_ptrlsts_.end();++iter)
    {
        if( obj_num_ < iter->size() ) obj_num_ = iter->size();
    }
}

void JRCSBox::update_color_label()
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
        Cube::colorByLabel((uint32_t*)vl.memptr(),vl.size(),label);
    }
}

void JRCSBox::init_from_boxes()
{
    std::vector<Cube::PtrLst>::iterator iter;
    MatPtrLst::iterator vviter = vvs_ptrlst_.begin();
    CMatPtrLst::iterator vciter = vcs_ptrlst_.begin();
    obj_prob_.clear();
    for(iter=cube_ptrlsts_.begin();iter!=cube_ptrlsts_.end();++iter)
    {
        if( iter->size()==obj_num_ && obj_prob_.empty() )
        {
            obj_prob_ = obj_prob_from_boxes(*iter,*vviter);
        }
        ++vviter;
        if(vviter==vvs_ptrlst_.end())break;
        ++vciter;
        if(vciter==vcs_ptrlst_.end())break;
    }
}

arma::fvec JRCSBox::obj_prob_from_boxes(const Cube::PtrLst& cube,const MatPtr& vv)
{
    for( Cube::PtrLst::const_iterator iter = cube.cbegin() ; iter!=cube.cend() ; ++iter )
    {
        ;
    }
}

}
