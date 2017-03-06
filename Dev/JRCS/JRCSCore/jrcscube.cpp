#include "jrcscube.h"
namespace JRCS{

Cube::Cube():plate_centroids_(3,plate_num_for_cube_,arma::fill::zeros),corners_(3,8,arma::fill::zeros),R_(3,3,arma::fill::eye),t_(3,arma::fill::zeros)
{
    xv_.reset(new arma::fmat(3,20));
    xn_.reset(new arma::fmat(3,20));
    xc_.reset(new arma::Mat<uint8_t>(3,20));
    t_ = obj_pos_;
    scale_r_ = arma::linspace<arma::fvec>(0.1,1.2,10);
}

Cube::Cube(
    const arma::fmat& v,
    const arma::fmat& n,
    const arma::Mat<uint8_t>& c,
    const arma::fvec& pos
):plate_centroids_(3,plate_num_for_cube_,arma::fill::zeros),corners_(3,8,arma::fill::zeros),R_(3,3,arma::fill::eye),t_(3,arma::fill::zeros)
{
    xv_.reset(new arma::fmat((float*)v.memptr(),v.n_rows,v.n_cols,false,true));
    xn_.reset(new arma::fmat((float*)n.memptr(),n.n_rows,n.n_cols,false,true));
    xc_.reset(new arma::Mat<uint8_t>((uint8_t*)c.memptr(),c.n_rows,c.n_cols,false,true));

    updateV2Centroids();

    updateV2Corners();

    t_ = arma::fvec(3,arma::fill::zeros);
    size_ = arma::max(arma::abs(corners_),1);
    obj_pos_ = pos;
    scale_r_ = arma::linspace<arma::fvec>(0.1,1.2,10);
}

void Cube::updateV2Centroids(void)
{
    for(int i = 0; i < xv_->n_cols ;i += point_num_for_plate_)
    {
        plate_centroids_.col(i) = arma::mean(xv_->cols(i,i+point_num_for_plate_-1));
    }
}

void Cube::updateV2Corners(void)
{
    ;
}

void Cube::updateCorners2V(void)
{
    for( int i = 0 ; i < xv_->n_cols ; ++i )
    {
        ;
    }
}

void Cube::translate(
        const arma::fvec& t,
        Cube& result
        )
{
    result.t_ = t_ + t;
    *result.xv_ = *xv_;
    result.xv_->each_col() += t;
    updateV2Centroids();
    updateV2Corners();
    if(this!=&result)
    {
        result.R_ = R_;
        *result.xc_ = *xc_;
        *result.xn_ = *xn_;
        result.size_ = size_;
        result.weighted_corners_ = weighted_corners_;
        result.obj_pos_ = obj_pos_ + t;
    }
}

void Cube::transform(
        const arma::fmat& R,
        const arma::fvec& t,
        Cube& result
        )
{
    result.R_ = R*R_;
    result.t_ = R*t_ + t;
    *result.xv_ = R*(*xv_);
    result.xv_->each_col() += t;
    *result.xn_ = R*(*xn_);
    updateV2Centroids();
    updateV2Corners();
    if(this!=&result)
    {
        *result.xc_ = *xc_;
        result.size_ = size_;
        result.weighted_corners_ = weighted_corners_;
        result.obj_pos_ = R*obj_pos_ + t;
    }
}

void Cube::scale(
        const arma::fvec& s,
        Cube& result
        )
{
    result.size_ = size_ % s;
    result.corners_ = R_.i()*corners_;
    result.corners_.each_col() %= s;
    result.corners_ = R_*result.corners_;
    updateCorners2V();
    updateV2Centroids();
    if(this!=&result)
    {
        result.t_ = t_;
        result.R_ = R_;
        *result.xc_ = *xc_;
        *result.xn_ = *xn_;
        result.weighted_corners_ = weighted_corners_;
        result.obj_pos_ = obj_pos_;
    }
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

void Cube::get_weighted_corners(
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

void Cube::accumulate(
        const arma::fmat& v,
        const arma::fmat& n,
        const arma::Mat<uint8_t>& c,
        const arma::vec alpha
        )
{
    ;
}

void Cube::start_accumulate(const int r,const int c,const int s,const int num)
{
    ;
}

void Cube::accumulate(const Cube&, const int i)
{
    ;
}

void Cube::fit(void)
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
            cube_ptrlst_[j]->transform(R,t,*cube_t_ptrlst_[i][j]);
        }
    }

    if(verbose_>1)std::cerr<<"calculate alpha"<<std::endl;
    arma::mat&  alpha = *alpha_ptrlst_[i];
    for(int c = 0 ; c < alpha.n_cols ; ++c )
    {
        alpha.col(c) = cube_t_ptrlst_[i][c]->get_dist2(vv_);
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
    assert(cube_t_ptrlst_[i].size()==alpha.n_cols);
    for(int c = 0 ; c < alpha.n_cols ; ++c )
    {
        cube_t_ptrlst_[i][c]->get_weighted_corners(vv_,alpha.col(c));
        cube_t_ptrlst_[i][c]->get_weighted_color(vv_,vc_);
    }

    if(verbose_>1)std::cerr<<"#2 calculating RT for each object"<<std::endl;
    for(int o = 0 ; o < obj_num_ ; ++o )
    {
        arma::fmat R(rt[o].R,3,3,false,true);
        arma::fvec t(rt[o].t,3,false,true);
        updateRTforObj(
                obj_range_[2*o],
                obj_range_[2*o+1],
                alpha_colsum,
                R,
                t,
                cube_t_ptrlst_[i]
                );
    }
    if(verbose_>1)std::cerr<<"#3 done RT for each object"<<std::endl;
    arma::mat alpha_v2(alpha.n_rows,alpha.n_cols);
    for(int c=0;c<alpha_v2.n_cols;++c)
    {
        alpha_v2.col(c) = cube_t_ptrlst_[i][c]->get_dist2(vv_);
    }
    vvar_.row(i) = arma::sum(alpha_v2%alpha);
    if(verbose_>1)std::cerr<<"#4 done var for each object"<<std::endl;
    for(int c=0;c<alpha.n_cols;++c)
    {
        cube_t_ptrlst_[i][c]->accumulate(vv_,vn_,vc_,alpha.col(c));
    }
    if(verbose_>1)std::cerr<<"#5 done accumulation"<<std::endl;
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
    arma::fmat A;
    arma::fmat U,V;
    arma::fvec s;
    arma::fmat dR;
    arma::fvec dt;

    arma::fmat _v(3,Cube::point_num_for_cube_);
    arma::fmat objv(3,Cube::point_num_for_cube_);
    int vi = 0;
    for(int p=start;p<=end;++p)
    {
        _v.col(vi) = cube_ptrlst[p]->weighted_corners_;
        objv.col(vi) = cube_ptrlst[p]->corners_;
        ++vi;
    }
    arma::fmat cv = _v.each_col() - arma::mean(_v,1);
    if(!_v.is_finite())
    {
        std::cerr<<iter_count_<<":"<<start<<"->"<<end<<":!v.is_finite()"<<std::endl;
    }

    arma::frowvec square_lambdav = arma::conv_to<arma::frowvec>::from(xv_invvar_.cols(start,end));
    square_lambdav %= colsum.cols(start,end);
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
    for(int p=start;p<=end;++p)
    {
        cube_ptrlst[p]->transform(dR,dt,*cube_ptrlst[p]);
    }
    //updating R T
    R = dR*R;
    t = dR*t + dt;
}

void JRCSCube::step_2(void)
{
    if(verbose_>1)std::cerr<<"Updating Latent Model"<<std::endl;

    #pragma omp parallel for
    for( int i=0 ; i < cube_ptrlst_.size() ; ++i )
    {
        const int r = cube_t_ptrlst_[0][i]->param_.n_rows;
        const int c = cube_t_ptrlst_[0][i]->param_.n_cols;
        const int s = cube_t_ptrlst_[0][i]->param_.n_slices;
        cube_ptrlst_[i]->start_accumulate(r,c,s,cube_t_ptrlst_.size());
        for(int j=0;j<cube_t_ptrlst_.size();++j)
        {
            cube_ptrlst_[i]->accumulate(*cube_t_ptrlst_[j][i],j);
        }
        cube_ptrlst_[i]->fit();
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

bool JRCSCube::isEnd_cube(void)
{
    if(iter_count_>=max_init_iter_)return true;
    else return false;
}

}
