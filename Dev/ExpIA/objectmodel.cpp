#include "objectmodel.h"
#include <nanoflann.hpp>
#include "common.h"
using namespace nanoflann;
ObjModel::ObjModel():
    GeoM_(new MeshBundle<DefaultMesh>),
    GeoLayout_(new MeshBundle<DefaultMesh>),
    ColorM_(new arma::gmm_diag())
{
    ;
}

void ObjModel::update(MeshBundle<DefaultMesh>::Ptr input)
{
    if(0==GeoM_->mesh_.n_vertices()){
        //init on first update
        GeoM_->mesh_.request_vertex_colors();
        GeoM_->mesh_.request_vertex_normals();
        for (DefaultMesh::VertexIter v_it = input->mesh_.vertices_begin();
             v_it != input->mesh_.vertices_end(); ++v_it)
        {
            DefaultMesh::VertexHandle new_v = GeoM_->mesh_.add_vertex(input->mesh_.point(*v_it));
            GeoM_->mesh_.set_color( new_v ,input->mesh_.color(*v_it));
        }
        GeoP_.resize(input->mesh_.n_vertices(),0.0);
        ColorP_.resize(input->mesh_.n_vertices(),0.0);
        return;
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
    arma::uvec indices(8);
    arma::fvec dists(8);
    arma::frowvec color_dists(8);
    arma::fmat neighbor_color;
    arma::fvec current_color;
    float* pdata = (float*)input->mesh_.points();
    uint8_t* cdata = (uint8_t*)input->mesh_.vertex_colors();
    arma::Mat<uint8_t> cmat(cdata,3,input->mesh_.n_vertices(),false,true);
    for(size_t index=0;index<input->mesh_.n_vertices();++index)
    {
        kdtree.knnSearch(&(pdata[3*index]),8,indices.memptr(),dists.memptr());
        space_density_(index) = dists(1);
        ColorArray::RGB2Lab(cmat.cols(indices),neighbor_color);
        ColorArray::RGB2Lab(cmat.col(index),current_color);
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
    MeshKDTreeInterface<DefaultMesh> pts(GeoM_->mesh_);
    KDTreeSingleIndexAdaptor<
            L2_Simple_Adaptor<float,MeshKDTreeInterface<DefaultMesh>>,
            MeshKDTreeInterface<DefaultMesh>,
            3,arma::uword>
            mkdtree(3,pts,KDTreeSingleIndexAdaptorParams(3));
    mkdtree.buildIndex();
    uint8_t* mcdata = (uint8_t*)GeoM_->mesh_.vertex_colors();
    arma::Mat<uint8_t> mcmat(mcdata,3,GeoM_->mesh_.n_vertices(),false,true);
    arma::fvec mc;
    arma::fvec c;
    size_t index = 0;
    for (DefaultMesh::VertexIter v_it = input->mesh_.vertices_begin();
         v_it != input->mesh_.vertices_end(); ++v_it)
    {
        mkdtree.knnSearch(&(pdata[3*index]),1,indices.memptr(),dists.memptr());
        float space_dist = dists(0);
        ColorArray::RGB2Lab(cmat.col(index),c);
        ColorArray::RGB2Lab(mcmat.col(indices(0)),mc);
        float color_dist = arma::norm( c - mc );
        if( (color_dist <= color_density_(index)) && (space_dist <= space_density_(index)) )
        {
            float wc = 1.0 / ( color_dist + 1.0 );
            float ws = 1.0 / ( space_dist + 1.0 );
            GeoP_[indices(0)] += wc;
            ColorP_[indices(0)] += ws;
        }else{
            DefaultMesh::VertexHandle new_v = GeoM_->mesh_.add_vertex(input->mesh_.point(*v_it));
            GeoM_->mesh_.set_color( new_v ,input->mesh_.color(*v_it));
            GeoP_.push_back(0.0);
            ColorP_.push_back(0.0);
        }
        ++index;
    }

}

void ObjModel::computeLayout()
{
    if(GeoP_.size()!=GeoM_->mesh_.n_vertices())std::logic_error("GeoP_.size()!=GeoM_->mesh_.n_vertices()");
    if(ColorP_.size()!=GeoM_->mesh_.n_vertices())std::logic_error("ColorP_.size()!=GeoM_->mesh_.n_vertices()");
    buildBB(GeoLayout_->mesh_);
    arma::fmat box((float*)GeoLayout_->mesh_.points(),3,8,false,true);
    arma::fmat pts((float*)GeoM_->mesh_.points(),3,GeoM_->mesh_.n_vertices(),false,true);
    arma::uvec indices;
    arma::fvec gp(GeoP_);
    arma::fvec cp(ColorP_);
    float gm = arma::max(gp);
    float cm = arma::max(cp);
    indices = arma::find( ( gp >= ( 0.6*gm ) ) && ( cp >= ( 0.6*cm ) ) );
    arma::fmat input = pts.cols(indices);
    arma::fmat R;
    arma::fvec t;
    std::vector<ObjModel::T::Ptr>::iterator iter;
    for(iter=GeoT_.begin();iter!=GeoT_.end();++iter)
    {
        if(0!=iter->use_count()&&(*iter))
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
    m = GeoM_->mesh_;
    arma::fmat V((float*)m.points(),3,m.n_vertices(),false,true);
    V = R*V + t;
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
    arma::fvec ColorP(ColorP_);
    if(!ColorP.save(path+"\\ColorP.fvec.arma",arma::arma_binary))return false;
    arma::fvec GeoP(GeoP_);
    if(!GeoP.save(path+"\\GeoP.fvec.arma",arma::arma_binary))return false;
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
    arma::fvec ColorP;
    if(!ColorP.load(path+"\\ColorP.fvec.arma"))return false;
    ColorP_ = arma::conv_to<std::vector<float>>::from(ColorP);
    arma::fvec GeoP;
    if(!GeoP.load(path+"\\GeoP.fvec.arma"))return false;
    GeoP_ = arma::conv_to<std::vector<float>>::from(GeoP);
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
