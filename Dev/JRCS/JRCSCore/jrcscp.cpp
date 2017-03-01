#include "jrcscp.h"
namespace JRCS {
CPlate::CPlate():Plate()
{
    ;
}

CPlate::CPlate(
    const arma::fmat& v,
    const arma::fmat& n,
    const arma::Mat<uint8_t>& c,
    const arma::fvec& pos
):Plate(v,n,c,pos)
{
    if(size_(2)==0.0)
    {
        type_ = H;
    }else{
        type_ = V;
    }
}

void CPlate::local_translate(
    const arma::fvec& t,
    Plate& result
)
{
    switch (type_) {
    case H:
        Plate::local_translate(t,result);
        break;
    case V:
        local_translate_v(t,result);
        break;
    }
}

void CPlate::scale(
    const arma::fvec &s,
    Plate &result
)
{
    switch (type_) {
    case H:
        Plate::scale(s,result);
        break;
    case V:
        scale_v(s,result);
        break;
    }
}

void CPlate::local_translate_v(const arma::fvec &t, Plate &result)
{
    arma::fvec dt = t;
    dt(2) = 0;
    *result.xv_ = *xv_;
    result.xv_->each_col() += dt;
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

void CPlate::scale_v(const arma::fvec& s, Plate& result)
{
    result.size_ = size_ % s;
    result.corners_ = R_.i()*corners_;
    result.corners_.each_col() %= s;
    result.corners_ = R_*result.corners_;
    *result.xv_ = result.corners_;
    result.xv_ -> each_col() += centroid_;
    //move the scaled plate so that it starts from floor
    //that is min(xv_->row(2)) == 0
    float z = arma::min(result.xv_->row(2));
    arma::fvec dt(3,arma::fill::zeros);
    dt(2) = obj_pos_(2) - z;
    result.xv_ -> each_col() += dt;
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

JRCSCP::JRCSCP():JRCSPrimitive()
{
    ;
}

void JRCSCP::initx(
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
                new CPlate(
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

void JRCSCP::prepare_primitive()
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
                        new CPlate(
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

}
