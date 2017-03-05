#include "jrcscube.h"
namespace JRCS{
Cube::Cube()
{
    ;
}

Cube::Cube(
    const arma::fmat& v,
    const arma::fmat& n,
    const arma::Mat<uint8_t>& c,
    const arma::fvec& pos
)
{
    ;
}

void Cube::translate(
        const arma::fvec& t,
        Cube& result
        )
{
    ;
}

void Cube::transform(
        const arma::fmat& R,
        const arma::fvec& t,
        Cube& result
        )
{
    ;
}

void Cube::scale(
        const arma::fvec& s,
        Cube& result
        )
{
    ;
}

arma::vec Cube::get_dist2(const arma::fmat& v)
{
    ;
}

arma::vec Cube::dist(
        const arma::fmat&,
        const arma::fvec&,
        arma::uword dim
        )
{
    ;
}

void Cube::get_weighted_centroid(
        const arma::fmat& v,
        const arma::vec &alpha
        )
{
    ;
}

void Cube::get_weighted_color(
        const arma::fmat& v,
        const arma::Mat<uint8_t>& c
        )
{
    ;
}

JRCSCube::JRCSCube():JRCSBilateral()
{

}

void JRCSCube::initx(
        const MatPtr& xv,
        const MatPtr& xn,
        const CMatPtr& xc
        )
{
    if(verbose_)std::cerr<<"JRCSCube start init x"<<std::endl;

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
    cube_ptrlst_.resize(obj_num_);
    arma::uword* pxr = (arma::uword*)obj_range_.data();
    arma::uword* pxr_s = pxr;
    *pxr = 0;
    for(int obj_idx = 0 ; obj_idx < obj_num_ ; ++ obj_idx )
    {
        int obj_size = Cube::point_num_for_cube_;
        if( pxr > pxr_s )*pxr = (*(pxr - 1))+1;
        *(pxr+1) = *pxr + 1 - 1;

        arma::fmat objv(pxv,3,obj_size,false,true);
        arma::fmat objn(pxn,3,obj_size,false,true);
        arma::Mat<uint8_t> objc(pxc,3,obj_size,false,true);

        arma::fvec pos = obj_pos_.col(obj_idx);
        reset_obj_vn(1.0,pos,objv,objn);
        reset_obj_c(objc);

        cube_ptrlst_[obj_idx].reset(
                    new Cube(
                        arma::fmat(pxv,3,Cube::point_num_for_cube_,false,true),
                        arma::fmat(pxn,3,Cube::point_num_for_cube_,false,true),
                        arma::Mat<uint8_t>(pxc,3,Cube::point_num_for_cube_,false,true),
                        obj_pos_.col(obj_idx)
                        )
                    );

        pxv += 3*Cube::point_num_for_cube_;
        pxn += 3*Cube::point_num_for_cube_;
        pxc += 3*Cube::point_num_for_cube_;

        pxr += 2;
        r_k -= obj_size;
    }
    if(verbose_)std::cerr<<"Done initx_"<<std::endl;
}

void JRCSCube::reset_obj_vn(
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

void JRCSCube::compute(void)
{
    if(verbose_)std::cerr<<"preparing"<<std::endl;
    prepare_cube();
    while(!isEnd_cube())
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
        finish_cube();
//        QThread::currentThread()->sleep(60);
    }
}

void JRCSCube::reset_alpha_cube()
{
    if(verbose_>0)std::cerr<<"allocating alpha"<<std::endl;
    int N = cube_ptrlst_.size();
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

void JRCSCube::reset_prob_cube()
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
    xv_invvar_ = arma::rowvec(cube_ptrlst_.size());
    xv_invvar_.fill(1.0/maxvar);

    vvar_ = arma::mat(vvs_ptrlst_.size(),cube_ptrlst_.size());

    x_p_ = arma::rowvec(cube_ptrlst_.size());
    x_p_.fill(1.0/float(cube_ptrlst_.size()));

    if(verbose_>0)std::cerr<<"done probability"<<std::endl;
}

void JRCSCube::prepare_cube()
{
    iter_count_ = 0;
    reset_alpha_cube();
    cube_t_ptrlst_.resize(vvs_ptrlst_.size());
    #pragma omp parallel for
    for( int i = 0 ; i < vvs_ptrlst_.size() ; ++i )
    {
        cube_t_ptrlst_[i].resize(cube_ptrlst_.size());
        float* pwv = (float*)wvs_ptrlst_[i]->memptr();
        float* pwn = (float*)wns_ptrlst_[i]->memptr();
        uint8_t* pwc = (uint8_t*)wcs_ptrlst_[i]->memptr();
        for(int j=0 ; j < cube_t_ptrlst_[i].size() ; ++j )
        {
            cube_t_ptrlst_[i][j].reset(
                        new Cube(
                            arma::fmat(pwv,3,Cube::point_num_for_cube_,false,true),
                            arma::fmat(pwn,3,Cube::point_num_for_cube_,false,true),
                            arma::Mat<uint8_t>(pwc,3,Cube::point_num_for_cube_,false,true),
                            arma::fvec(3,arma::fill::zeros)
                            )
                        );
            pwv += 3*Cube::point_num_for_cube_;
            pwn += 3*Cube::point_num_for_cube_;
            pwc += 3*Cube::point_num_for_cube_;
        }
    }
    reset_prob_cube();
}

void JRCSCube::finish_cube()
{
    ++iter_count_;
}

void JRCSCube::step_1(int i)
{
    ;
}

void JRCSCube::updateRTforObj(
        const int start,
        const int end,
        arma::frowvec& colsum,
        arma::fmat& R,
        arma::fvec& t,
        Cube::PtrLst cube_ptrlst
        )
{
    ;
}

void JRCSCube::step_2(void)
{
    ;
}

bool JRCSCube::isEnd_cube(void)
{
    if(iter_count_>=max_init_iter_)return true;
    else return false;
}

}
