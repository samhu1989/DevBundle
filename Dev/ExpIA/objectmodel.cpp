#include "objectmodel.h"
#include <nanoflann.hpp>
#include "common.h"
#include "featurecore.h"
using namespace nanoflann;
ObjModel::ObjModel():
    GeoM_(new MeshBundle<DefaultMesh>),
    GeoLayout_(new MeshBundle<DefaultMesh>)
{
    ;
}

void ObjModel::init(arma::fmat&X)
{
    if(0==GeoM_->mesh_.n_vertices()){
        //init on first update
        accu_color_ = arma::mat(3,X.n_cols,arma::fill::zeros);
        accu_count_ = 0;
        DistP_  = arma::fvec(X.n_cols,arma::fill::zeros);
        ColorP_ = arma::fvec(X.n_cols,arma::fill::zeros);
        NormP_ = arma::fvec(X.n_cols,arma::fill::zeros);
        GeoM_->mesh_.request_vertex_colors();
        GeoM_->mesh_.request_vertex_normals();
        for( size_t i = 0 ; i < X.n_cols ; ++ i )
        {
            GeoM_->mesh_.add_vertex(DefaultMesh::Point(X(0,i),X(1,i),X(2,i)));
        }
        Feature::computePointNormal<DefaultMesh>(GeoM_->mesh_,0.0,10);
        return;
    }
}

void ObjModel::updateColor(MeshBundle<DefaultMesh>::Ptr input)
{
    MeshKDTreeInterface<DefaultMesh> points(input->mesh_);
    KDTreeSingleIndexAdaptor<
            L2_Simple_Adaptor<float,MeshKDTreeInterface<DefaultMesh>>,
            MeshKDTreeInterface<DefaultMesh>,
            3,arma::uword>
            kdtree(3,points,KDTreeSingleIndexAdaptorParams(2));
    kdtree.buildIndex();
    arma::uvec indices(1);
    arma::fvec dists(1);
    arma::fvec current_color;
    arma::fvec current_normal;
    float* pdata = (float*)GeoM_->mesh_.points();
    uint8_t* cdata = (uint8_t*)input->mesh_.vertex_colors();
    arma::Mat<uint8_t> cmat(cdata,3,input->mesh_.n_vertices(),false,true);
    for(size_t index=0;index<GeoM_->mesh_.n_vertices();++index)
    {
        kdtree.knnSearch(&(pdata[3*index]),1,indices.memptr(),dists.memptr());
        double w = 1.0 / ( 1.0 + dists(0) );
        current_color = arma::conv_to<arma::fvec>::from(cmat.col(indices(0)));
        accu_color_.col(index) += arma::conv_to<arma::vec>::from(w*current_color);
        DistP_[index] += w;
    }
    ++ accu_count_;
}

void ObjModel::finishColor()
{
    #pragma omp for
    for(size_t index=0;index<GeoM_->mesh_.n_vertices();++index)
    {
        accu_color_.col(index) /= DistP_(index);
    }
    arma::Mat<uint8_t> mc((uint8_t*)GeoM_->mesh_.vertex_colors(),3,GeoM_->mesh_.n_vertices(),false,true);
    mc = arma::conv_to<arma::Mat<uint8_t>>::from(accu_color_);
    #pragma omp for
    for(size_t index=0;index<GeoM_->mesh_.n_vertices();++index)
    {
        DistP_(index) /= accu_count_;
    }
    accu_count_ = 0;
}

void ObjModel::updateWeight(MeshBundle<DefaultMesh>::Ptr input)
{
    MeshKDTreeInterface<DefaultMesh> points(input->mesh_);
    KDTreeSingleIndexAdaptor<
            L2_Simple_Adaptor<float,MeshKDTreeInterface<DefaultMesh>>,
            MeshKDTreeInterface<DefaultMesh>,
            3,arma::uword>
            kdtree(3,points,KDTreeSingleIndexAdaptorParams(2));
    kdtree.buildIndex();
    arma::uvec indices(1);
    arma::fvec dists(1);
    arma::fvec current_color;
    arma::fvec neighbor_color;
    arma::fvec current_normal;
    arma::fvec neighbor_normal;

    uint8_t* mcdata = (uint8_t*)GeoM_->mesh_.vertex_colors();
    arma::Mat<uint8_t> mcmat(mcdata,3,GeoM_->mesh_.n_vertices(),false,true);
    float* mndata = (float*)GeoM_->mesh_.vertex_normals();
    arma::fmat mnmat(mndata,3,GeoM_->mesh_.n_vertices(),false,true);

    float* pdata = (float*)GeoM_->mesh_.points();
    uint8_t* cdata = (uint8_t*)input->mesh_.vertex_colors();
    arma::Mat<uint8_t> cmat(cdata,3,input->mesh_.n_vertices(),false,true);
    float* ndata = (float*)input->mesh_.vertex_normals();
    arma::fmat nmat(ndata,3,input->mesh_.n_vertices(),false,true);

    for(size_t index=0;index<GeoM_->mesh_.n_vertices();++index)
    {
        kdtree.knnSearch(&(pdata[3*index]),1,indices.memptr(),dists.memptr());
        ColorArray::RGB2Lab(cmat.col(indices(0)),current_color);
        ColorArray::RGB2Lab(mcmat.col(index),neighbor_color);
        current_normal = nmat.col(indices(0));
        neighbor_normal = mnmat.col(index);
        double color_dist = arma::norm(current_color - neighbor_color);
        double cw = 1.0 / ( 1.0 + color_dist );
        ColorP_(index) += cw;
        if(current_normal.is_finite()&&neighbor_normal.is_finite())
        {
            double nw = std::abs(arma::dot(current_normal,neighbor_normal));
            NormP_(index) += nw;
        }
    }
    ++ accu_count_;
}

void ObjModel::finishWeight()
{
    #pragma omp for
    for(size_t index=0;index<GeoM_->mesh_.n_vertices();++index)
    {
        ColorP_(index) /= accu_count_;
        NormP_(index) /= accu_count_;
    }
    accu_count_ = 0;
}

void ObjModel::computeLayout()
{
    if(DistP_.size()!=GeoM_->mesh_.n_vertices())std::logic_error("GeoP_.size()!=GeoM_->mesh_.n_vertices()");
    if(NormP_.size()!=GeoM_->mesh_.n_vertices())std::logic_error("NormP_.size()!=GeoM_->mesh_.n_vertices()");
    if(ColorP_.size()!=GeoM_->mesh_.n_vertices())std::logic_error("ColorP_.size()!=GeoM_->mesh_.n_vertices()");
    buildBB(GeoLayout_->mesh_);
    arma::fmat box((float*)GeoLayout_->mesh_.points(),3,8,false,true);
    arma::fmat pts((float*)GeoM_->mesh_.points(),3,GeoM_->mesh_.n_vertices(),false,true);
    arma::uvec indices;
    float dm = arma::median(DistP_);
    float cm = arma::median(ColorP_);
    float nm = arma::median(NormP_);
    indices = arma::find( ( DistP_ >= dm ) && ( ColorP_ >= cm ) && ( NormP_ >= nm ) );
    if(indices.is_empty())
    {
       indices = arma::find(  ( ColorP_ >= cm ) );
    }
    arma::fmat input = pts.cols(indices);
    arma::fmat R;
    arma::fvec t;
    std::vector<ObjModel::T::Ptr>::iterator iter;
    for(iter=GeoT_.begin();iter!=GeoT_.end();++iter)
    {
        if((*iter)&&0!=iter->use_count())
        {
            R = arma::fmat((*iter)->R,3,3,false,true);
            t = arma::fvec((*iter)->t,3,false,true);
            break;
        }
    }
    input.each_col() -= t;
    input = arma::inv(R)*input;
    get3DMBB(input,2,box);
    box = R*box ;
    box.each_col() += t;
}

bool ObjModel::transform(DefaultMesh& m,uint32_t T_index)
{
    ObjModel::T::Ptr& T_ptr = GeoT_[T_index];
    if(T_ptr&&0!=T_ptr.use_count())
    {
        arma::fmat R(T_ptr->R,3,3,false,true);
        arma::fvec t(T_ptr->t,3,false,true);
        return transform(m,R,t);
    }else return false;
}

bool ObjModel::transform(DefaultMesh& m,arma::fmat& R,arma::fvec& t)
{
    m.request_vertex_normals();
    m.request_vertex_colors();
    m = GeoM_->mesh_;
    arma::fmat V((float*)m.points(),3,m.n_vertices(),false,true);
    V.each_col() -= t;
    V = arma::inv(R)*V;
    arma::fmat N((float*)m.vertex_normals(),3,m.n_vertices(),false,true);
    N = arma::inv(R)*N;
    return true;
}

bool ObjModel::save(const std::string& path)
{
    OpenMesh::IO::Options opt;
    opt+=OpenMesh::IO::Options::Binary;
    opt+=OpenMesh::IO::Options::VertexColor;
    opt+=OpenMesh::IO::Options::VertexNormal;
    GeoM_->mesh_.request_vertex_normals();
    if(!OpenMesh::IO::write_mesh(GeoM_->mesh_,path+"\\GeoM.ply",opt,10)){
        std::cerr<<"can't save to:"<<path+"\\GeoM.ply"<<std::endl;
        return false;
    }
    opt-=OpenMesh::IO::Options::VertexColor;
    opt-=OpenMesh::IO::Options::VertexNormal;
    if(!OpenMesh::IO::write_mesh(GeoLayout_->mesh_,path+"\\GeoLayout.ply",opt,10)){
        std::cerr<<"can't save to:"<<path+"\\GeoLayout.ply"<<std::endl;
        return false;
    }
    if(!ColorP_.save(path+"\\ColorP.fvec.arma",arma::arma_binary))return false;
    if(!NormP_.save(path+"\\NormP.fvec.arma",arma::arma_binary))return false;
    if(!DistP_.save(path+"\\DistP.fvec.arma",arma::arma_binary))return false;
    //saving Ts;
    arma::uvec T_indices(GeoT_.size(),arma::fill::zeros);
    std::vector<T::Ptr>::iterator iter;
    size_t index = 0;
    size_t T_cnt = 0;
    for(iter=GeoT_.begin();iter!=GeoT_.end();++iter)
    {
        T::Ptr& ptr = *iter;
        if(ptr&&0!=ptr.use_count())
        {
            T_indices(index) = 1;
            ++ T_cnt;
        }
        ++index;
    }
    if(!T_indices.save(path+"\\T_indices.uvec.arma",arma::arma_binary))return false;
    arma::fcube T_cube(3,4,T_cnt);
    index = 0;
    for(iter=GeoT_.begin();iter!=GeoT_.end();++iter)
    {
        T::Ptr& ptr = *iter;
        if(ptr&&0!=ptr.use_count())
        {
            arma::fmat R(ptr->R,3,3,false,true);
            arma::fmat t(ptr->t,3,1,false,true);
            T_cube.slice(index) = arma::join_rows(R,t);
            ++ index;
        }
    }
    if(!T_cube.save(path+"\\T_cube.fcube.arma",arma::arma_binary))return false;
    return true;
}

bool ObjModel::load(const std::string& path)
{
    OpenMesh::IO::Options opt;
    opt+=OpenMesh::IO::Options::Binary;
    opt+=OpenMesh::IO::Options::VertexColor;
    opt+=OpenMesh::IO::Options::VertexNormal;
    GeoM_->mesh_.request_vertex_normals();
    GeoM_->mesh_.request_vertex_colors();
    if(!OpenMesh::IO::read_mesh(GeoM_->mesh_,path+"\\GeoM.ply",opt))return false;
    GeoLayout_->mesh_.request_vertex_normals();
    GeoLayout_->mesh_.request_vertex_colors();
    if(!OpenMesh::IO::read_mesh(GeoLayout_->mesh_,path+"\\GeoLayout.ply",opt))return false;
    if(!ColorP_.load(path+"\\ColorP.fvec.arma"))return false;
    if(!ColorP_.is_finite())
    {
        std::cerr<<"Infinite in ColorP"<<std::endl;
    }
    if(!NormP_.load(path+"\\NormP.fvec.arma"))return false;
    if(!NormP_.is_finite())
    {
        std::cerr<<"Infinite in NormP"<<std::endl;
    }
    if(!DistP_.load(path+"\\DistP.fvec.arma"))return false;
    if(!DistP_.is_finite())
    {
        std::cerr<<"Infinite in DistP"<<std::endl;
    }
    //loading Ts;
    arma::uvec T_indices;
    if(!T_indices.load(path+"\\T_indices.uvec.arma"))return false;
    arma::fcube T_cube;
    if(!T_cube.load(path+"\\T_cube.fcube.arma"))return false;
    if( T_cube.n_cols != 4 || T_cube.n_rows != 3 )return false;
    GeoT_.resize(T_indices.size());
    arma::uvec indices = arma::find(T_indices==1);
    arma::uword index = 0;
    arma::uvec::iterator iter;
    for(iter=indices.begin();iter!=indices.end();++iter)
    {
        T::Ptr& ptr = GeoT_[*iter];
        ptr = std::make_shared<T>();
        arma::fmat R(ptr->R,3,3,false,true);
        arma::fvec t(ptr->t,3,false,true);
        R = T_cube.slice(index).cols(0,2);
        t = T_cube.slice(index).col(3);
        ++index;
    }
    return true;
}
