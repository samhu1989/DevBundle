#include "objectmodel.h"
#include <nanoflann.hpp>
#include "common.h"
#include "featurecore.h"
#include "filter.h"
using namespace nanoflann;
ObjModel::ObjModel():
    GeoM_(new MeshBundle<DefaultMesh>),
    GeoLayout_(new MeshBundle<DefaultMesh>)
{
    ;
}

void ObjModel::init(arma::fmat&X)
{
    initX_ = X;
    GeoM_->mesh_.request_vertex_colors();
    GeoM_->mesh_.request_vertex_normals();
}

void ObjModel::updateFullModel(MeshBundle<DefaultMesh>::Ptr input)
{
    if(0==FullM_.n_vertices())
    {
        FullM_ = input->mesh_;
    }
    //evaluate local density of current input
    arma::fvec space_density_(input->mesh_.n_vertices());
    arma::fvec color_density_(input->mesh_.n_vertices());
    MeshKDTreeInterface<DefaultMesh> points(input->mesh_);
    KDTreeSingleIndexAdaptor<
            L2_Simple_Adaptor<float,MeshKDTreeInterface<DefaultMesh>>,
            MeshKDTreeInterface<DefaultMesh>,
            3,arma::uword>
            kdtree(3,points,KDTreeSingleIndexAdaptorParams(9));
    kdtree.buildIndex();
    size_t k = std::min(size_t(8),input->mesh_.n_vertices());
    arma::uvec indices(k);
    arma::fvec dists(k);
    arma::frowvec color_dists;
    arma::fmat neighbor_color;
    arma::fvec current_color;
    float* pdata = (float*)input->mesh_.points();
    uint8_t* cdata = (uint8_t*)input->mesh_.vertex_colors();
    arma::Mat<uint8_t> cmat(cdata,3,input->mesh_.n_vertices(),false,true);
    for(size_t index=0;index<input->mesh_.n_vertices();++index)
    {
        kdtree.knnSearch(&(pdata[3*index]),k,indices.memptr(),dists.memptr());
        space_density_(index) = dists(1);
        neighbor_color = arma::conv_to<arma::fmat>::from(cmat.cols(indices));
        current_color = arma::conv_to<arma::fvec>::from(cmat.col(index));
        neighbor_color.each_col() -= current_color;
        color_dists = arma::sum(arma::square(neighbor_color));
        arma::frowvec tmp = arma::sort(color_dists);
        color_density_(index) = tmp(1);
    }
    //merge a point to current model
    //if it is within the local density range;
    //and add one to confidence( also counts )
    //add a point to current model
    //if it is not within the local density range;
    MeshKDTreeInterface<DefaultMesh> pts(FullM_);
    KDTreeSingleIndexAdaptor<
            L2_Simple_Adaptor<float,MeshKDTreeInterface<DefaultMesh>>,
            MeshKDTreeInterface<DefaultMesh>,
            3,arma::uword>
            mkdtree(3,pts,KDTreeSingleIndexAdaptorParams(3));
    mkdtree.buildIndex();
    size_t index = 0;
    for (DefaultMesh::VertexIter v_it = input->mesh_.vertices_begin();
         v_it != input->mesh_.vertices_end(); ++v_it)
    {
        mkdtree.knnSearch(&(pdata[3*index]),1,indices.memptr(),dists.memptr());
        float space_dist = dists(0);
        if( space_dist >= space_density_(index) )
        {
            DefaultMesh::VertexHandle new_v = FullM_.add_vertex(input->mesh_.point(*v_it));
            FullM_.set_color( new_v ,input->mesh_.color(*v_it));
        }
        ++index;
    }
}

void ObjModel::updateModel(MeshBundle<DefaultMesh>::Ptr input)
{
    MeshKDTreeInterface<DefaultMesh> points(input->mesh_);
    KDTreeSingleIndexAdaptor<
            L2_Simple_Adaptor<float,MeshKDTreeInterface<DefaultMesh>>,
            MeshKDTreeInterface<DefaultMesh>,
            3,arma::uword>
            kdtree(3,points,KDTreeSingleIndexAdaptorParams(3));
    kdtree.buildIndex();

    size_t N = initX_.n_cols;
    float* pdata = (float*)initX_.memptr();
    float* idata = (float*)input->mesh_.points();
    size_t k = std::min(size_t(6),input->mesh_.n_vertices());
    arma::uvec indices(k);
    arma::fvec dists(k);
    //init on first update
    for( size_t index = 0 ; index < N ; ++ index )
    {
        kdtree.knnSearch(&(pdata[3*index]),k,indices.memptr(),dists.memptr());
        for(size_t ii=0;ii<k;++ii)
        {
            GeoM_->mesh_.add_vertex(
                    DefaultMesh::Point(
                        idata[3*indices(ii)],
                        idata[3*indices(ii)+1],
                        idata[3*indices(ii)+2])
                );
        }
    }
}

void ObjModel::finishModel()
{
    Filter::OctreeGrid<DefaultMesh> filter;
    filter.set_seed_resolution(0.02);
    filter.extract(GeoM_->mesh_);
    Feature::computePointNormal<DefaultMesh>(GeoM_->mesh_,0.0,10);
    accu_color_ = arma::mat(3,GeoM_->mesh_.n_vertices(),arma::fill::zeros);
    accu_count_ = 0;
    DistP_  = arma::fvec(GeoM_->mesh_.n_vertices(),arma::fill::zeros);
    ColorP_ = arma::fvec(GeoM_->mesh_.n_vertices(),arma::fill::zeros);
    NormP_ = arma::fvec(GeoM_->mesh_.n_vertices(),arma::fill::zeros);
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
        current_color = arma::conv_to<arma::fvec>::from(cmat.col(indices(0)));
        neighbor_color = arma::conv_to<arma::fvec>::from(mcmat.col(index));
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
    float dm = arma::mean(DistP_);
    float cm = arma::mean(ColorP_);
    float nm = arma::mean(NormP_);
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

void ObjModel::computeFullLayout()
{
    buildBB(FullLayout_);
    arma::fmat box((float*)FullLayout_.points(),3,8,false,true);
    arma::fmat pts((float*)FullM_.points(),3,FullM_.n_vertices(),true,true);
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
    pts.each_col() -= t;
    pts = arma::inv(R)*pts;
    get3DMBB(pts,2,box);
    box = R*box ;
    box.each_col() += t;
}

bool ObjModel::fullLayout(arma::fmat& object_layout,int32_t t_index)
{
    object_layout = arma::fmat((float*)FullLayout_.points(),3,FullLayout_.n_vertices(),true,true);
    ObjModel::T::Ptr t_ptr;
    if(t_index>=0&&t_index<GeoT_.size())
    {
        t_ptr = GeoT_[t_index];
        if(t_ptr&&0!=t_ptr.use_count())
        {
            arma::fmat R(t_ptr->R,3,3,false,true);
            arma::fvec t(t_ptr->t,3,false,true);
            object_layout.each_col() -= t;
            object_layout = arma::inv(R)*object_layout;
        }else return false;
    }
    return true;
}

bool ObjModel::fullModel(DefaultMesh& object_mesh,int32_t t_index)
{
    object_mesh = FullM_;
    arma::fmat model_mat((float*)object_mesh.points(),3,object_mesh.n_vertices(),false,true);
    ObjModel::T::Ptr t_ptr;
    if(t_index>=0&&t_index<GeoT_.size())
    {
        t_ptr = GeoT_[t_index];
        if(t_ptr&&0!=t_ptr.use_count())
        {
            arma::fmat R(t_ptr->R,3,3,false,true);
            arma::fvec t(t_ptr->t,3,false,true);
            model_mat.each_col() -= t;
            model_mat = arma::inv(R)*model_mat;
        }else return false;
    }
    return true;
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
    if(!OpenMesh::IO::write_mesh(FullM_,path+"\\FullM.ply",opt,10)){
        std::cerr<<"can't save to:"<<path+"\\FullM.ply"<<std::endl;
        return false;
    }
    opt-=OpenMesh::IO::Options::VertexColor;
    opt-=OpenMesh::IO::Options::VertexNormal;
    if(!OpenMesh::IO::write_mesh(GeoLayout_->mesh_,path+"\\GeoLayout.ply",opt,10)){
        std::cerr<<"can't save to:"<<path+"\\GeoLayout.ply"<<std::endl;
        return false;
    }
    if(!OpenMesh::IO::write_mesh(FullLayout_,path+"\\FullLayout.ply",opt,10)){
        std::cerr<<"can't save to:"<<path+"\\FullLayout.ply"<<std::endl;
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
    FullM_.request_vertex_normals();
    FullM_.request_vertex_colors();
    if(!OpenMesh::IO::read_mesh(FullM_,path+"\\FullM.ply",opt))return false;
    GeoLayout_->mesh_.request_vertex_normals();
    GeoLayout_->mesh_.request_vertex_colors();
    if(!OpenMesh::IO::read_mesh(GeoLayout_->mesh_,path+"\\GeoLayout.ply",opt))return false;
    FullLayout_.request_vertex_normals();
    FullLayout_.request_vertex_colors();
    if(!OpenMesh::IO::read_mesh(FullLayout_,path+"\\FullLayout.ply",opt))return false;
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
