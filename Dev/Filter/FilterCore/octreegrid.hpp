#include "octreegrid.h"
namespace Filter {
template<typename M>
bool OctreeGrid<M>::initCompute()
{
    if(0==seed_resolution_)return false;
    if(!std::isfinite(seed_resolution_))return false;
    return true;
}
template<typename M>
void OctreeGrid<M>::applyFilter(M&mesh)
{
    M result;
    if(mesh.has_vertex_colors())
    {
        result.request_vertex_colors();
    }
    OctreeMesh octree_m(mesh);
    arma::fmat p_mat((float*)mesh.points(),3,mesh.n_vertices(),false,true);
    arma::Mat<uint8_t> c_mat((uint8_t*)mesh.vertex_colors(),3,mesh.n_vertices(),false,true);
    Octree tree;
    tree.initialize(octree_m,unibn::OctreeParams(float(seed_resolution_)));
    size_t N = tree.getLeafNum();
    for( size_t index = 0 ; index < N ; ++ index )
    {
        arma::uvec pIndices;
        if(0==tree.getPointofLeafAt(index,pIndices))continue;
        arma::fvec p = arma::mean( p_mat.cols(pIndices) , 1 );
        arma::fmat cs = arma::conv_to<arma::fmat>::from( c_mat.cols(pIndices) );
        arma::Col<uint8_t> c = arma::conv_to<arma::Col<uint8_t>>::from( arma::mean(cs,1) );
        MeshVertexHandle new_v = result.add_vertex(MeshPoint(p(0),p(1),p(2)));
        result.set_color( new_v , MeshColor(c(0),c(1),c(2)) );
    }
    mesh = result;
}
}
