#include "cube.h"
#include <QColor>
namespace Common {
std::vector<arma::uvec> Cube::c4v_;
arma::fvec Cube::scale_r_;
const uint32_t Cube::point_num_for_plate_ = 4;
const uint32_t Cube::plate_num_for_cube_ = 5;
const uint32_t Cube::point_num_for_cube_ = 20;
__gnu_cxx::hash_map<uint32_t,uint32_t,std::hash<uint32_t>,Cube::color_equal_to> Cube::color_label_;
__gnu_cxx::hash_map<uint32_t,uint32_t> Cube::label_color_;

uint32_t Cube::colorFromLabel(uint32_t label)
{
    if(label==0)return QColor("white").rgba();
    uint32_t c;
    if( Cube::label_color_.end() == Cube::label_color_.find(label) )
    {
        c = ColorArray::rand_color();
        while( Cube::color_label_.end() != Cube::color_label_.find(c) )//if color is duplicated rand another one
        {
            c = ColorArray::rand_color();
        }
        Cube::label_color_[label] = c;
        Cube::color_label_[c] = label;
    }else{
        c = Cube::label_color_[label];
    }
    return c;
}

void Cube::colorByLabel(uint32_t label)
{
    QColor color(colorFromLabel(label));
    arma::Col<uint8_t> x = {color.red(),color.green(),color.blue()};
    xc_->each_col() = x;
}

void Cube::colorByLabel(uint32_t* c,arma::uword size,arma::uvec& label)
{
    if( label.size() < size ){
        std::cerr<<"label.size() < size"<<std::endl;
        return;
    }
    for(arma::uvec::iterator iter=label.begin() ; iter != label.end() ; ++iter )
    {
        *c = colorFromLabel(*iter);
        ++c;
    }
}

Cube::PtrLst Cube::newCubes(DefaultMesh& m, uint32_t N)
{
    Cube::PtrLst r(N);
    uint32_t mN = 20;

    for(int n=0 ; n < N ; ++n)
    {
        for(int i=0 ; i < 5 ; ++i)
        {
            std::vector<DefaultMesh::VertexHandle>  face_vhandles_a,face_vhandles_b;
            face_vhandles_a.push_back(m.add_vertex(DefaultMesh::Point(0,0,0)));
            face_vhandles_a.push_back(m.add_vertex(DefaultMesh::Point(0,0,0)));
            face_vhandles_a.push_back(m.add_vertex(DefaultMesh::Point(0,0,0)));
            face_vhandles_b.push_back(m.add_vertex(DefaultMesh::Point(0,0,0)));
            face_vhandles_b.push_back(face_vhandles_a[0]);
            face_vhandles_b.push_back(face_vhandles_a[2]);

            m.add_face(face_vhandles_a);
            m.add_face(face_vhandles_b);
        }
    }

    m.request_face_normals();
    m.request_face_colors();
    m.request_vertex_normals();
    m.request_vertex_colors();

    float* pp = (float*)m.points();
    float* pn = (float*)m.vertex_normals();
    uint8_t* pc = (uint8_t*)m.vertex_colors();

    arma::fvec pos = {0,0,0};
    for(int n = 0  ;n < N ; ++n )
    {
        arma::fmat xv(pp,3,mN,false,true);
        arma::fmat xn(pn,3,mN,false,true);
        arma::Mat<uint8_t> xc(pc,3,mN,false,true);

        xv = {
            //0  1  2  3  4  5  6  7  8  9  10 11 12 13 14 15 16 17 18 19
            {-1, 1, 1,-1, 1,-1,-1, 1, 1, 1, 1, 1,-1,-1,-1,-1,-1,-1, 1, 1},
            { 1, 1, 1, 1,-1,-1,-1,-1, 1,-1,-1, 1,-1, 1, 1,-1, 1,-1,-1, 1},
            { 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 1, 1}
        };

        xv *= 0.5;

        xn = {
            { 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1,-1,-1,-1,-1, 0, 0, 0, 0},
            { 1, 1, 1, 1,-1,-1,-1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1}
        };

        arma::Col<uint8_t> xctmp = {137,157,192};
        xc.each_col() = xctmp;
        r[n].reset( new Cube(xv,xn,xc,pos) );

        pp += 3*mN;
        pn += 3*mN;
        pc += 3*mN;
    }
    return r;
}

Cube::Ptr Cube::newCube(void)
{
    arma::fmat xv(3,20);
    arma::fmat xn(3,20);
    arma::Mat<uint8_t> xc(3,20);
    xv = {
        //0  1  2  3  4  5  6  7  8  9  10 11 12 13 14 15 16 17 18 19
        {-1, 1, 1,-1, 1,-1,-1, 1, 1, 1, 1, 1,-1,-1,-1,-1,-1,-1, 1, 1},
        { 1, 1, 1, 1,-1,-1,-1,-1, 1,-1,-1, 1,-1, 1, 1,-1, 1,-1,-1, 1},
        { 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 1, 1}
    };
    xv *= 0.5;
    xn = {
        { 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1,-1,-1,-1,-1, 0, 0, 0, 0},
        { 1, 1, 1, 1,-1,-1,-1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1}
    };
    arma::Col<uint8_t> xctmp = {137,157,192};
    xc.each_col() = xctmp;
    arma::fvec pos = {0,0,0};
    return Cube::Ptr(new Cube(xv,xn,xc,pos));
}

Cube::Ptr Cube::newCube(DefaultMesh& mesh)
{
    uint32_t N = mesh.n_vertices();

    for(int i=0 ; i < 5 ; ++i)
    {
        std::vector<DefaultMesh::VertexHandle>  face_vhandles_a,face_vhandles_b;
        face_vhandles_a.push_back(mesh.add_vertex(DefaultMesh::Point(0,0,0)));
        face_vhandles_a.push_back(mesh.add_vertex(DefaultMesh::Point(0,0,0)));
        face_vhandles_a.push_back(mesh.add_vertex(DefaultMesh::Point(0,0,0)));
        face_vhandles_b.push_back(mesh.add_vertex(DefaultMesh::Point(0,0,0)));
        face_vhandles_b.push_back(face_vhandles_a[0]);
        face_vhandles_b.push_back(face_vhandles_a[2]);

        mesh.add_face(face_vhandles_a);
        mesh.add_face(face_vhandles_b);
    }

    mesh.request_face_normals();
    mesh.request_face_colors();
    mesh.request_vertex_normals();
    mesh.request_vertex_colors();

    float* p = (float*)mesh.points();
    float* n = (float*)mesh.vertex_normals();
    uint8_t* c = (uint8_t*)mesh.vertex_colors();

    p += 3*N;
    n += 3*N;
    c += 3*N;

    arma::fmat xv(p,3,20,false,true);
    arma::fmat xn(n,3,20,false,true);
    arma::Mat<uint8_t> xc(c,3,20,false,true);

    xv = {
        //0  1  2  3  4  5  6  7  8  9  10 11 12 13 14 15 16 17 18 19
        {-1, 1, 1,-1, 1,-1,-1, 1, 1, 1, 1, 1,-1,-1,-1,-1,-1,-1, 1, 1},
        { 1, 1, 1, 1,-1,-1,-1,-1, 1,-1,-1, 1,-1, 1, 1,-1, 1,-1,-1, 1},
        { 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 1, 1}
    };

    xv *= 0.5;

    xn = {
        { 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1,-1,-1,-1,-1, 0, 0, 0, 0},
        { 1, 1, 1, 1,-1,-1,-1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1}
    };

    arma::Col<uint8_t> xctmp = {137,157,192};
    xc.each_col() = xctmp;

    arma::fvec pos = {0,0,0};

    return Cube::Ptr(new Cube(xv,xn,xc,pos));
}

Cube::Cube():plate_centroids_(3,plate_num_for_cube_,arma::fill::zeros),corners_(3,8,arma::fill::zeros),R_(3,3,arma::fill::eye),t_(3,arma::fill::zeros)
{

    xv_.reset(new arma::fmat(3,20));
    xn_.reset(new arma::fmat(3,20));
    xc_.reset(new arma::Mat<uint8_t>(3,20));
    t_ = obj_pos_;
    if(scale_r_.empty())scale_r_ = arma::linspace<arma::fvec>(0.5,1.5,11) ;
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

    if(scale_r_.empty())scale_r_ = arma::linspace<arma::fvec>(0.5,1.5,11) ;
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

arma::fvec Cube::bottom_pos()
{
    arma::fvec result(3,arma::fill::zeros);
    result.rows(0,1) = arma::mean(corners_.rows(0,1),1);
    result.row(2) = arma::min(corners_.row(2),1);
    return result;
}

void Cube::translate(
        const arma::fvec& t,
        Cube& result
        )
{
    result.t_ = t_ + t;
    *result.xv_ = *xv_;
    result.xv_->each_col() += t;
    result.obj_pos_ = obj_pos_ + t;
    result.updateV2Centroids();
    result.updateV2Corners();
    if(this!=&result)
    {
        result.R_ = R_;
        *result.xc_ = *xc_;
        *result.xn_ = *xn_;
        result.size_ = size_;
        result.plate_zero_dim_ = plate_zero_dim_;
        result.weighted_corners_ = weighted_corners_;
    }
}

void Cube::rotate(
        const arma::fmat& R,
        Cube& result
        )
{
    arma::fvec t = t_;
    translate(-t,result);
    transform(R,t,result);
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
    result.obj_pos_ = R*obj_pos_ + t;

    result.updateV2Centroids();
    result.updateV2Corners();

    if(this!=&result)
    {
        *result.xc_ = *xc_;
        result.size_ = size_;
        result.plate_zero_dim_ = plate_zero_dim_;
        result.weighted_corners_ = weighted_corners_;
    }
}

void Cube::updateScale()
{
    arma::fmat c = corners_.each_col() - obj_pos_;
    c = R_.i()*c ;
    arma::fmat wc = weighted_corners_.each_col() - obj_pos_;
    wc = R_.i()*wc ;

    arma::fmat A = wc*arma::pinv(c);
    arma::fmat U,V;
    arma::fvec s;
    arma::svd(U,s,V,A,"std");
    arma::fvec ss = s;
    arma::fmat absU = arma::abs(U);
    arma::fvec dim0 = absU.col(0);
    arma::fvec dim1 = absU.col(1);
    arma::fvec dim2 = absU.col(2);

    arma::uword maxi;
    dim0.max(maxi);
    ss(maxi) = s(0);
    dim1.max(maxi);
    ss(maxi) = s(1);
    dim2.max(maxi);
    ss(maxi) = s(2);

    ss(0) = size_(0)*ss(0) > 0.01 ? ss(0) : 1.0;
    ss(1) = size_(1)*ss(1) > 0.01 ? ss(1) : 1.0;
    ss(2) = size_(2)*ss(2) > 0.01 ? ss(2) : 1.0;

//    std::cerr<<"ss:"<<ss<<std::endl;

    scale(ss,*this);
}

void Cube::scaleTo(
        const arma::fvec& news
        )
{
    arma::fvec s = news / size_;
    scale(s,*this);
}

void Cube::scaleTo(
        const arma::fvec& news,
        Cube& result
        )
{
    arma::fvec s = news / size_;
    scale(s,result);
}

void Cube::scale(
        const arma::fvec& s,
        Cube& result
        )
{
    result.size_ = size_ % s;
    result.corners_ =  corners_.each_col() - obj_pos_;
    result.corners_ = R_.i()*result.corners_;
    result.corners_.each_col() %= s;
    result.corners_ = R_*result.corners_;
    result.corners_.each_col() += obj_pos_;

    result.updateCorners2V();
    result.updateV2Centroids();

    if(this!=&result)
    {
        result.t_ = t_;
        result.R_ = R_;
        *result.xc_ = *xc_;
        *result.xn_ = *xn_;
        result.plate_zero_dim_ = plate_zero_dim_;
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

arma::vec Cube::get_dist2_box(const arma::fmat& v)
{
    arma::vec dist(v.n_cols,arma::fill::zeros);
    arma::uvec out_idx = outside(v);
    dist(out_idx) = get_dist2(v.cols(out_idx));
}

arma::uvec Cube::outside(const arma::fmat& v)
{
    std::vector<arma::uword> idx;
    idx.reserve(v.n_cols);
    arma::fvec center = arma::mean(corners_,1);
    center = R_.i()*( center - t_ );


    return arma::uvec(idx);
}

arma::uvec Cube::inside(const arma::fmat& v)
{
    ;
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
    mean();
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
        scale_size(0) = size_(0)*scale_r_(i) > 0.01 ? scale_r_(i) : 1.0;
        scale_size(1) = size_(1)*scale_r_(j) > 0.01 ? scale_r_(j) : 1.0;
        scale_size(2) = size_(2)*scale_r_(k) > 0.01 ? scale_r_(k) : 1.0;
        std::cerr<<"Cube::fit scale:"<<scale_size.t()<<std::endl;
        scale(scale_size,*this);
    }
    param_.clear();
}



}
