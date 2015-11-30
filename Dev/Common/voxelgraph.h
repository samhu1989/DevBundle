#ifndef VOXELGRAPH_H
#define VOXELGRAPH_H
#include <armadillo>
template <typename M>
class VoxelGraph
{
public:
    typedef M Mesh;
    VoxelGraph(M&m):Ref_(m){}
    bool save(const std::string&);
    bool load(const std::string&);
    void sv2pix(arma::uvec& sv,arma::uvec& pix);//supervoxel label to pixel label
    void match(M&mesh,std::vector<float>& gscore, std::vector<float>& cscore,arma::vec&score);
    double voxel_similarity(size_t v1,size_t v2);

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
