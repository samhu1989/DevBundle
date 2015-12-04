#ifndef MESHTYPE_H
#define MESHTYPE_H
#include <memory>
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>
#include <OpenMesh/Tools/Utils/StripifierT.hh>
#include <armadillo>
#include "MeshColor.h"
#include "voxelgraph.h"
struct Traits : public OpenMesh::DefaultTraits
{
  HalfedgeAttributes(OpenMesh::Attributes::PrevHalfedge);
};


typedef OpenMesh::PolyMesh_ArrayKernelT<Traits>  DefaultMesh;

template <typename M>
class MeshBundle
{
public:
    typedef typename std::shared_ptr<MeshBundle<M>> Ptr;
    typedef typename std::vector<Ptr> PtrList;
    typedef typename std::shared_ptr<M> MeshPtr;
    MeshBundle():
        custom_color_(mesh_),
        graph_(mesh_),
        strips_(mesh_)
    {}
    MeshPtr mesh_ptr(){return std::shared_ptr<M>(&mesh_);}
    std::string name_;
    M                           mesh_;
    MeshColor<M>        custom_color_;
    VoxelGraph<M>              graph_;
    OpenMesh::StripifierT<M>  strips_;
};

template <typename IM,typename OM>
inline void extractMesh(const IM& input,arma::uvec& indices,OM& output)
{
    typedef typename OM::Point P;
    arma::fmat pts((float*)input.points(),3,input.n_vertices(),false,true);
    for( int i=0 ; i < indices.size() ; ++i )
    {
        arma::uword index = indices(i);
        output.add_vertex(P(pts(0,index),pts(1,index),pts(2,index)));
    }
    output.request_vertex_normals();
    output.request_vertex_colors();
    arma::Mat<uint8_t> cIn((uint8_t*)input.vertex_colors(),3,input.n_vertices(),false,true);
    arma::Mat<uint8_t> cOut((uint8_t*)output.vertex_colors(),3,output.n_vertices(),false,true);
    cOut = cIn.cols(indices);
}



#endif // MESHTYPE

