#include "voxelgraph.h"
#include "KDtree.hpp"
#include "nanoflann.hpp"
template <typename M>
bool VoxelGraph<M>::save(const std::string&path)
{
    if(!voxel_centers.save(path+"\\centers.fvec.arma",arma::arma_binary))return false;
    if(!voxel_neighbors.save(path+"\\neighbors.Mat_uint16_t.arma",arma::arma_binary))return false;
    if(!voxel_label.save(path+"\\labels.uvec.arma",arma::arma_binary))return false;
    return true;
}
template <typename M>
bool VoxelGraph<M>::load(const std::string&path)
{
    if(!voxel_centers.load(path+"\\centers.fvec.arma"))return false;
    if(!voxel_neighbors.load(path+"\\neighbors.Mat_uint16_t.arma"))return false;
    if(!voxel_label.load(path+"\\labels.uvec.arma"))return false;
    return true;
}
template <typename M>
void VoxelGraph<M>::sv2pix(arma::uvec& sv,arma::uvec& pix)
{
    if(sv.size()!=voxel_centers.n_cols){
        std::cerr<<"Can't translate a supervoxel label that was not based on this graph"<<std::endl;
        std::logic_error("sv.size()!=voxel_centers.n_cols");
    }
    if(pix.size()!=voxel_label.size())pix = arma::uvec(voxel_label.size(),arma::fill::zeros);
    else pix.fill(0);
    for(int l = 1 ; l <= voxel_centers.n_cols ; ++l )
    {
        arma::uvec indices = arma::find(voxel_label==l);
        pix(indices).fill( sv(l-1) );
    }
}
using namespace nanoflann;
template <typename M>
void VoxelGraph<M>::match(
        M&mesh,
        std::vector<float>&gscore,
        std::vector<float>&cscore,
        std::vector<float>&score
        )
{
    MeshKDTreeInterface<M> points(mesh);
    KDTreeSingleIndexAdaptor<
            L2_Simple_Adaptor<float,MeshKDTreeInterface<M>>,
            MeshKDTreeInterface<M>,
            3,arma::uword>
            kdtree(3,points,KDTreeSingleIndexAdaptorParams(9));
    kdtree.buildIndex();
    SearchParams param;
    size_t sv_N = voxel_centers.n_cols;
    arma::uvec search_idx(5);
    arma::fvec search_dist(5);
    std::vector<float>sv_match_score(sv_N,0.0f);
    std::vector<float>sv_geo_score(sv_N,0.0f);
    std::vector<float>sv_color_score(sv_N,0.0f);
    score.resize(sv_N);
    float *pts = (float*)Ref_.points();
    for( size_t p_i = 0 ; p_i < voxel_centers.n_cols ; ++ p_i )
    {
        kdtree.knnSearch(&pts[3*p_i],5,search_idx.memptr(),search_dist.memptr());
    }
    for(size_t sv_i = 0 ; sv_i < voxel_centers.n_cols ; ++ sv_i )
    {
        score[sv_i] = sv_match_score[sv_i]*sv_geo_score[sv_i]*sv_color_score[sv_i];
    }
}

template <typename M>
double VoxelGraph<M>::voxel_similarity(size_t v1,size_t v2)
{
    ;
}
