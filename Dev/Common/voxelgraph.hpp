#include "voxelgraph.h"
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
