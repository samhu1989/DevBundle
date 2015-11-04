#ifndef SUPERVOXELCLUSTERING_H
#define SUPERVOXELCLUSTERING_H
#include <memory>
#include <functional>
#include <armadillo>
#include <ext/hash_map>
#include "segmentationbase.h"
#include "common.h"
namespace Segmentation{
template<typename M>
class SuperVoxel
{
public:
    typedef std::shared_ptr<SuperVoxel> Ptr;
    SuperVoxel(const M&mesh):
        mesh_(mesh),
        voxels_((float*)mesh.points(),3,mesh.n_vertices(),false,true),
        normals_((float*)mesh.vertex_normals(),3,mesh.n_vertices(),false,true),
        colors_((uint8_t*)mesh.vertex_colors(),3,mesh.n_vertices(),false,true)
    {}
    inline arma::fvec& centroid(){return centeroid_;}
    inline arma::fvec& normal(){return normal_;}
    inline arma::fvec& color() {return rgb_;}
    const arma::fmat& voxels(){return voxels_.cols(indices_);}
    const arma::fmat& normals(){return normals_.cols(indices_);}
    const arma::Mat<uint8_t>& colors(){return colors_.cols(indices_);}
    arma::uvec& indices(){return indices_;}
    const arma::uvec& indices()const{return indices_;}
protected:
    arma::fvec centeroid_;
    arma::fvec normal_;
    arma::fvec rgb_;
    const arma::fmat normals_;
    const arma::fmat voxels_;
    const arma::Mat<uint8_t> colors_;
    arma::uvec indices_;
private:
    const M& mesh_;
};

template<typename M>
class DefaultVoxelDistFunctor
{
public:
    float dist(const typename SuperVoxel<M>::Ptr& v1,const typename SuperVoxel<M>::Ptr& v2,float seed_res)
    {
        float spatial_dist = arma::norm( v1->centroid() - v2->centroid() ) / seed_res;
        float color_dist = arma::norm( v1->color() - v2->color() ) / 255.0f;
        float cos_angle_normal = 1.0f - std::abs ( arma::dot( v1->normal() , v2->normal() ) );
        return  cos_angle_normal * normal_importance_ + color_dist * color_importance_+ spatial_dist * spatial_importance_;
    }
    float normal_importance_;
    float color_importance_;
    float spatial_importance_;
};

template<typename M>
class SuperVoxelClustering:public SegmentationBase
{
public:
    typedef std::function<float(const typename SuperVoxel<M>::Ptr&,const typename SuperVoxel<M>::Ptr&)> DistFunc;
    typedef __gnu_cxx::hash_map<uint64_t,typename SuperVoxel<M>::Ptr> SuperVoxels;
    typedef __gnu_cxx::hash_multimap<uint64_t,uint64_t> SuperVoxelAdjacency;
    typedef std::shared_ptr<M> MeshPtr;
    SuperVoxelClustering(float v_res,float seed_res);
    SuperVoxelClustering(float v_res,float seed_res,bool);
    virtual ~SuperVoxelClustering();

    void setDistFunctor(DistFunc&);
    void input(MeshPtr mesh);
    void extract(arma::uvec&labels);
    void extract(SuperVoxels&supervoxel_clusters);
    void getSupervoxelAdjacency(SuperVoxelAdjacency&label_adjacency)const;
protected:
    bool initCompute();
    void deinitCompute();
    void selectInitialSupervoxelSeeds(arma::uvec&indices);
    void createSupervoxels(arma::uvec& seed_indices);
    void expandSupervoxels(int depth);
    void makeSupervoxels(SuperVoxels&);
protected:
    MeshPtr input_;
private:
    float seed_resolution_;
    float voxel_resolution_;
    DistFunc vDist_;
};
}
#include "supervoxelclustering.hpp"
#endif // SUPERVOXELCLUSTERING_H
