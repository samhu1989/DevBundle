#include "jrcscube.h"
#include <QThread>
namespace JRCS{

std::vector<arma::uvec> Cube::c4v_;

Cube::Cube():plate_centroids_(3,plate_num_for_cube_,arma::fill::zeros),corners_(3,8,arma::fill::zeros),R_(3,3,arma::fill::eye),t_(3,arma::fill::zeros)
{

    xv_.reset(new arma::fmat(3,20));
    xn_.reset(new arma::fmat(3,20));
    xc_.reset(new arma::Mat<uint8_t>(3,20));
    t_ = obj_pos_;
    scale_r_ = arma::linspace<arma::fvec>(0.5,1.5,11) ;
}

Cube::Cube(
    const arma::fmat& v,
    const arma::fmat& n,
    const arma::Mat<uint8_t>& c,
    const arma::fvec& pos
):plate_centroids_(3,plate_num_for_cube_,arma::fill::zeros),corners_(3,8,arma::fill::zeros),R_(3,3,arma::fill::eye),t_(3,arma::fill::zeros)
{
    assert(plate_num_for_cube_*point_num_for_plate_==point_num_for_cube_);
    xv_.reset(new arma::fmat((float*)v.memptr(),v.n_rows,v.n_cols,false,true));
    xn_.reset(new arma::fmat((float*)n.memptr(),n.n_rows,n.n_cols,false,true));
    xc_.reset(new arma::Mat<uint8_t>((uint8_t*)c.memptr(),c.n_rows,c.n_cols,false,true));

    updateV2Centroids();
    updateV2Corners();
    updateZeroDim();

    t_ = arma::fvec(3,arma::fill::zeros);
    updateCorners2Size();
    obj_pos_ = pos;
    scale_r_ = arma::linspace<arma::fvec>(0.5,1.5,11);
}

void Cube::updateCorners2Size(void)
{
    arma::fmat tmp = corners_;
    tmp.each_col() -= arma::mean(tmp,1);
    size_ = arma::max( arma::abs(tmp) , 1 );
}

void Cube::updateZeroDim(void)
{
    plate_zero_dim_ = arma::uvec( plate_num_for_cube_ );
    for(int i = 0; i < plate_num_for_cube_ ; ++i )
    {
        int start = point_num_for_plate_*i;
        int end = point_num_for_plate_*(i+1) - 1;
        arma::fmat tmp = xv_->cols(start,end);
        arma::fvec tmpc = arma::mean(tmp,1);
        tmp.each_col() -= tmpc;
        arma::fvec maxm = arma::max(arma::abs(tmp),1);
        for(arma::uword idx=0;idx<3;++idx)
        {
            if( maxm(idx) < std::numeric_limits<float>::epsilon() )
            {
                plate_zero_dim_(i) = idx;
                break;
            }
        }
    }
}

void Cube::updateV2Centroids(void)
{
    for(int i = 0; i < plate_num_for_cube_ ; ++i )
    {
        int start = point_num_for_plate_*i;
        int end = point_num_for_plate_*(i+1) - 1;
        plate_centroids_.col(i) = arma::mean(xv_->cols(start,end),1);
    }
}

void Cube::updateV2Corners(void)
{
    if(Cube::c4v_.empty())
    {
        Cube::c4v_.resize(8);
        Cube::c4v_[0] = {0,13,16};
        Cube::c4v_[1] = {1,8,19};
        Cube::c4v_[2] = {2,11};
        Cube::c4v_[3] = {3,14};
        Cube::c4v_[4] = {5,12,17};
        Cube::c4v_[5] = {4,9,18};
        Cube::c4v_[6] = {7,10};
        Cube::c4v_[7] = {6,15};
    }
    uint32_t i = 0;
    for(std::vector<arma::uvec>::iterator iter = c4v_.begin() ; iter != c4v_.end()  ; ++iter )
    {
        corners_.col(i) = arma::mean(xv_->cols(*iter),1);
        ++i;
    }
}

void Cube::updateCorners2V(void)
{
    if(Cube::c4v_.empty())
    {
        Cube::c4v_.resize(8);
        Cube::c4v_[0] = {0,13,16};
        Cube::c4v_[1] = {1,8,19};
        Cube::c4v_[2] = {2,11};
        Cube::c4v_[3] = {3,14};
        Cube::c4v_[4] = {5,12,17};
        Cube::c4v_[5] = {4,9,18};
        Cube::c4v_[6] = {7,10};
        Cube::c4v_[7] = {6,15};
    }
    uint32_t i = 0;
    for(std::vector<arma::uvec>::iterator iter = c4v_.begin() ; iter != c4v_.end()  ; ++iter )
    {
        for(int c=0;c<iter->size();++c)
        {
            xv_->col((*iter)(c)) = corners_.col(i);
        }
        ++i;
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
    result.updateV2Centroids();
    result.updateV2Corners();
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

    result.updateV2Centroids();
    result.updateV2Corners();

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

    result.updateCorners2V();
    result.updateV2Centroids();

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
    arma::mat dists(v.n_cols,plate_centroids_.n_cols);
    for(uint32_t i = 0 ; i < plate_centroids_.n_cols ; ++i )
    {
        dists.col(i) = get_dist2_for_plate(v,plate_centroids_.col(i),plate_zero_dim_(i));
    }
    return arma::min(dists,1);
}

arma::vec Cube::get_dist2_for_plate(
        const arma::fmat& v,
        const arma::fvec& c,
        arma::uword zero_dim
        )
{
    arma::fmat tv = v.each_col() - t_;
    arma::fmat invR = R_.i();
    arma::fvec o = c;
    o -= t_;
    assert(invR.is_finite());
    tv = invR*tv;
    o = invR*o;
    return dist(tv,o,zero_dim,0)+dist(tv,o,zero_dim,1)+dist(tv,o,zero_dim,2);
}

arma::vec Cube::dist(
        const arma::fmat& v,
        const arma::fvec& origin,
        arma::uword zero_dim,
        arma::uword dim
        )
{
    arma::vec result(v.n_cols,arma::fill::zeros);
    arma::frowvec vdim = v.row(dim);
    arma::uvec idx_dim;
    if(dim==zero_dim)
    {
        idx_dim = arma::find( arma::abs( vdim - origin(dim) ) > 0 );
        result(idx_dim) = arma::square( arma::abs( arma::conv_to<arma::vec>::from( vdim.cols(idx_dim)) - origin(dim) ) );
    }else{
        idx_dim = arma::find( arma::abs( vdim - origin(dim) ) > size_(dim) );
        result(idx_dim) = arma::square( arma::abs( arma::conv_to<arma::vec>::from( vdim.cols(idx_dim)) - origin(dim) ) - size_(dim) );
    }
    return result;
}

void Cube::get_weighted_corners(
        const arma::fmat& v,
        const arma::vec &alpha
        )
{
    weighted_corners_ = corners_;
    for( uint32_t i = 0 ; i < corners_.n_cols  ; ++i )
    {
        arma::fvec c = corners_.col(i);
        arma::rowvec a = arma::trunc_exp( arma::conv_to<arma::rowvec>::from( - arma::sum(  arma::square( v.each_col() - c ) ) ) );
        a %= alpha.t();
        a += std::numeric_limits<double>::epsilon();
        weighted_corners_.col(i) = arma::conv_to<arma::fvec>::from( v*( a.t() ) / arma::accu(a) );
    }
}

void Cube::get_weighted_color(
        const arma::fmat& v,
        const arma::Mat<uint8_t>& c
        )
{
    arma::fmat w(v.n_cols,xv_->n_cols);
    for( int i = 0 ; i < w.n_cols; ++i )
    {
        w.col(i) = arma::trunc_exp(-1.0*arma::sum(arma::square( v.each_col() - xv_->col(i) ))).t();
    }
    w += std::numeric_limits<float>::epsilon();
    arma::frowvec sum = arma::sum(w);
    w.each_row() /= sum;
    arma::fmat fc = arma::conv_to<arma::fmat>::from(c)*w;
    *xc_ = arma::conv_to<arma::Mat<uint8_t>>::from(fc);
}

void Cube::accumulate(
        const arma::fmat& v,
        const arma::fmat& n,
        const arma::Mat<uint8_t>& c,
        const arma::vec alpha
        )
{
    arma::fvec scale_size(3,arma::fill::zeros);
    if(param_.empty())param_ = arma::fcube(scale_r_.size(),scale_r_.size(),scale_r_.size(),arma::fill::zeros);
    else param_.fill(0.0);
    for(int i = 0;i<scale_r_.size();++i)
    {
        for(int j=0;j<scale_r_.size();++j)
        {
            for(int k=0;k<scale_r_.size();++k)
            {
                arma::fmat tmpv((float*)xv_->memptr(),xv_->n_rows,xv_->n_cols,true,true);
                arma::fmat tmpn((float*)xn_->memptr(),xn_->n_rows,xn_->n_cols,true,true);
                arma::Mat<uint8_t> tmpc((uint8_t*)xc_->memptr(),xc_->n_rows,xc_->n_cols,true,true);
                Cube::Ptr tmp_cube(new Cube(tmpv,tmpn,tmpc,obj_pos_));
                scale_size(0) = scale_r_(i);
                scale_size(1) = scale_r_(j);
                scale_size(2) = scale_r_(k);
                scale(scale_size,*tmp_cube);
                arma::vec dist2 = tmp_cube->get_dist2(v);
                dist2 = arma::trunc_exp( - dist2 ) % alpha ;
                param_(i,j,k) = arma::accu( dist2 );
            }
        }
    }
}

void Cube::start_accumulate(const int r,const int c,const int s,const int num)
{
    param_ = arma::fcube(r,c,s,arma::fill::zeros);
    param_vec_.reset(new arma::fvec(param_.memptr(),param_.size(),false,true));
    param_mat_.reset(new arma::fmat(param_.size(),num));
}

void Cube::accumulate(const Cube& p, const int i)
{
    param_mat_->col(i) = arma::fvec((float*)p.param_.memptr(),p.param_.size(),false,true);
}

void Cube::median(void)
{
    (*param_vec_) = arma::median(*param_mat_,1);
}

void Cube::mean(void)
{
    (*param_vec_) = arma::mean(*param_mat_,1);
}

void Cube::fit(void)
{
    median();
    arma::uword i=0,j=0,k=0;
    float fitscore = param_.max(i,j,k);
    std::cerr<<"fit score="<<fitscore<<"@("<<i<<","<<j<<","<<k<<")"<<std::endl;
    param_(i,j,k) = 0.0;
    arma::uword a=0,b=0,c=0;
    float sscore = param_.max(a,b,c);
    std::cerr<<"fit second score="<<sscore<<"@("<<a<<","<<b<<","<<c<<")"<<std::endl;
    if(fitscore<=std::numeric_limits<float>::epsilon())
    {
        std::cerr<<"maximum fit score below zero:"<<fitscore<<std::endl;
    }
    if( fitscore > 0.0 )
    {
        arma::fvec scale_size(3,arma::fill::zeros);
        scale_size(0) = scale_r_(i);
        scale_size(1) = scale_r_(j);
        scale_size(2) = scale_r_(k);
        std::cerr<<"Cube::fit scale:"<<scale_size.t()<<std::endl;
        scale(scale_size,*this);
    }
    param_.clear();
}

JRCSCube::JRCSCube():JRCSBilateral()
{
    ;
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
        reset_obj_vn(0.8,pos,objv,objn);
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

    if(verbose_)std::cerr<<"JRCSCube end init x"<<std::endl;
}

void JRCSCube::reset_obj_vn(
        float radius,
        arma::fvec& pos,
        arma::fmat& ov,
        arma::fmat& on
        )
{
//    if(verbose_)std::cerr<<"JRCSCube::reset_obj_vn"<<std::endl;
    arma::fmat v = {
        //0  1  2  3  4  5  6  7  8  9  10 11 12 13 14 15 16 17 18 19
        {-1, 1, 1,-1, 1,-1,-1, 1, 1, 1, 1, 1,-1,-1,-1,-1,-1,-1, 1, 1},
        { 1, 1, 1, 1,-1,-1,-1,-1, 1,-1,-1, 1,-1, 1, 1,-1, 1,-1,-1, 1},
        { 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 1, 1}
    };
    //for Cube
    Cube::c4v_.resize(8);
    Cube::c4v_[0] = {0,13,16};
    Cube::c4v_[1] = {1,8,19};
    Cube::c4v_[2] = {2,11};
    Cube::c4v_[3] = {3,14};
    Cube::c4v_[4] = {5,12,17};
    Cube::c4v_[5] = {4,9,18};
    Cube::c4v_[6] = {7,10};
    Cube::c4v_[7] = {6,15};
    ov = v;
    ov *= radius;
    ov.each_col() += pos;
    arma::fmat n = {
        { 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1,-1,-1,-1,-1, 0, 0, 0, 0},
        { 1, 1, 1, 1,-1,-1,-1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1}
    };
    on = n;
//    if(verbose_)std::cerr<<"JRCSCube::reset_obj_vn:end"<<std::endl;
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
            Cube& tc = *cube_t_ptrlst_[i][j];
            arma::vec v = tc.get_dist2(*wvs_ptrlst_[i]);
            ColorArray::colorfromValue((ColorArray::RGB888*)wcs_ptrlst_[i]->memptr(),wcs_ptrlst_[i]->n_cols,arma::sqrt(v));
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
    if(verbose_>1)std::cerr<<"#1 calculate weighted cube corners"<<std::endl;
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
    if(verbose_)std::cerr<<"updating RT"<<std::endl;
    arma::fmat A;
    arma::fmat U,V;
    arma::fvec s;
    arma::fmat dR;
    arma::fvec dt;

    arma::fmat _v(3,8);
    arma::fmat objv(3,8);

    assert(start==end);
    if(start!=end)
    {
        std::cerr<<"start should be equal to end for cube"<<std::endl;
    }

    _v = cube_ptrlst[start]->weighted_corners_;
    objv = cube_ptrlst[start]->corners_;

    arma::fmat cv = _v.each_col() - arma::mean(_v,1);
    if(!_v.is_finite())
    {
        std::cerr<<iter_count_<<":"<<start<<"->"<<end<<":!v.is_finite()"<<std::endl;
    }

    arma::fmat tmpv = objv.each_col() - arma::mean(objv,1) ;
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
            dt = arma::mean(tmp,1);
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
            dt = arma::mean(tmp,1);
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
    if(verbose_)std::cerr<<"done updating RT"<<std::endl;
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
