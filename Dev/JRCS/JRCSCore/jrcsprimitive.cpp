#include "jrcsprimitive.h"
namespace JRCS{

Plate::Plate()
{
    ;
}

Plate::Plate(
    const arma::fmat& v,
    const arma::fmat& n,
    const arma::Mat<uint8_t>& c
)
{
    xv_.reset(new arma::fmat(v.memptr(),v.n_rows,v.n_cols));
    xn_.reset(new arma::fmat(n.memptr(),v.n_rows,v.n_cols));
    xc_.reset(new arma::Mat<uint8_t>(c.memptr(),c.n_rows,c.n_cols));
}

void Plate::transform(
        const arma::fmat& R,
        const arma::fvec& t,
        Plate&
        )
{
    ;
}

arma::vec Plate::get_alpha(
        const arma::fmat& v
        )
{
    ;
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

    max_obj_radius_ = 1.0;// for randomly reset rt

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
        reset_obj_vn(0.5,pos,objv,objn);
        reset_obj_c(objc);

        for(int j=0 ; j < plate_num_for_obj_; ++j)
        {
            plate_ptrlst_[obj_idx*plate_num_for_obj_+j].reset(
                new Plate(
                    arma::fmat(pxv,3,point_num_for_plate_,false,true),
                    arma::fmat(pxn,3,point_num_for_plate_,false,true),
                    arma::Mat<uint8_t>(pxc,3,point_num_for_plate_,false,true)
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
        if(verbose_)std::cerr<<"step 1"<<std::endl;
        #pragma omp parallel for
        for( int i=0 ; i < vvs_ptrlst_.size() ; ++i )
        {
            step_1(i);
        }
        if(verbose_)std::cerr<<"step 2"<<std::endl;
        step_2();
        finish_primitive();
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
        for(int j=0 ; j < plate_t_ptrlst_[i].size() ; ++j )
        {
            plate_t_ptrlst_[i][j].reset(new Plate());
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
        alpha.col(c) = plate_t_ptrlst_[i][c]->get_alpha(vv_);
    }

    alpha = arma::trunc_exp(alpha);
    alpha.each_row() %=  arma::pow(xv_invvar_,1.5);

    if(verbose_>1)std::cerr<<"normalize alpha"<<std::endl;
    alpha += std::numeric_limits<double>::epsilon(); //add eps for numeric stability
    arma::vec alpha_rowsum = ( 1.0 + beta_ ) * arma::sum(alpha,1);
    alpha.each_col() /= alpha_rowsum;
}

void JRCSPrimitive::step_2(void)
{
    ;
}

bool JRCSPrimitive::isEnd_primitive(void)
{
    if(iter_count_>=1)return true;
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
