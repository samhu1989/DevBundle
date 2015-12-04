#ifndef VOXELGRAPH_H
#define VOXELGRAPH_H
#include <armadillo>
template <typename M>
class VoxelGraph
{
public:
    typedef M Mesh;
    VoxelGraph(M&m):Ref_(m){}
    bool empty(void)
    {
        return ( Ref_.n_vertices() != voxel_label.size() ) || ( 0==voxel_centers.n_cols );
    }

    bool save(const std::string&);
    bool load(const std::string&);
    void sv2pix(arma::uvec& sv,arma::uvec& pix);//supervoxel label to pixel label
    void sv2pix(arma::Col<uint32_t>&,arma::Col<uint32_t>&);//supervoxel color to pixel color
    void match(
            M&mesh,std::vector<float>& gscore,
            std::vector<float>& cscore,
            arma::vec&score,
            double dist_th=0.05,
            double color_var=30.0
            );
    double voxel_similarity(size_t v1,size_t v2,double dist_th=0.05,double color_var=30.0);
    inline void get_XYZLab(arma::fmat&voxels,const arma::uvec&indices = arma::uvec());
    inline void get_Lab(arma::fmat&voxels,const arma::uvec&indices = arma::uvec());

    arma::fmat voxel_centers;
    arma::Mat<uint8_t> voxel_colors;
    arma::uvec voxel_size;
    arma::uvec voxel_label;//start from one
    arma::Mat<uint16_t> voxel_neighbors;
private:
    const Mesh& Ref_;//the ref Mesh that bound with this custom color
};
#include "voxelgraph.hpp"
#endif // VOXELGRAPH_H