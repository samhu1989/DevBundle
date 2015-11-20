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
            kdtree(3,points,KDTreeSingleIndexAdaptorParams(4));
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
    indices = arma::find( ( gp >= 0.5*gm ) && ( cp >= 0.5*cm ) );
    arma::fmat input = pts.cols(indices);
    get3DMBB(input,2,box);
}
