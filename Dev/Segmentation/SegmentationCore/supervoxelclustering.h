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
class SuperVoxel;

template<typename M>
class SuperVoxelClustering;

template<typename M>
class Voxel
{
public:
    typedef std::shared_ptr<Voxel> Ptr;
    Voxel(const M&mesh,arma::uvec& indices):
        mesh_(mesh),
        points_((float*)mesh.points(),3,mesh.n_vertices(),false,true),
        normals_((float*)mesh.vertex_normals(),3,mesh.n_vertices(),false,true),
        colors_((uint8_t*)mesh.vertex_colors(),3,mesh.n_vertices(),false,true),
        indices_(indices)
    {
        xyz_ = arma::mean(points_,1);
        normal_ = arma::mean(normals_,1);
        arma::fmat colors = arma::conv_to<arma::fmat>::from(colors_);
        rgb_ = arma::mean(colors,1);
    }

    Voxel(const M& mesh):
        xyz_(arma::fvec(3,arma::fill::zeros)),
        rgb_(arma::fvec(3,arma::fill::zeros)),
        normal_(arma::fvec(3,arma::fill::zeros)),
        indices_(0),
        mesh_(mesh)
    {}

    inline arma::fvec& xyz(){return xyz_;}
    inline arma::fvec& normal(){return normal_;}
    inline arma::fvec& color() {return rgb_;}

    const arma::fmat& voxels()const{return points_.cols(indices_);}
    const arma::fmat& normals()const{return normals_.cols(indices_);}
    const arma::Mat<uint8_t>& colors()const{return colors_.cols(indices_);}
    const arma::uvec& indices()const{return indices_;}
    std::vector<Ptr> neighbors_;
    typename SuperVoxel<M>::Ptr parent_;
    float dist2parent_;
protected:

    arma::fvec xyz_;
    arma::fvec normal_;
    arma::fvec rgb_;

    const arma::fmat normals_;
    const arma::fmat points_;
    const arma::Mat<uint8_t> colors_;
    const arma::uvec& indices_;
private:
    const M& mesh_;
};

template<typename M>
class SuperVoxel
{
public:
    typedef typename std::vector<typename Voxel<M>::Ptr> VoxelVec;
    typedef typename std::vector<typename Voxel<M>::Ptr>::iterator VoxelIter;
    typedef std::shared_ptr<SuperVoxel> Ptr;
    typedef SuperVoxelClustering<M> Parent;
    SuperVoxel(const Parent& p):parent_(p){centroid_.reset(new Voxel<M>(*p.input_));}
    size_t size(){return voxels.size();}
    void updateCentroid();
    void expand();
    void removeVoxel(typename Voxel<M>::Ptr);
    void addVoxel(typename Voxel<M>::Ptr);
private:
    typename Voxel<M>::Ptr centroid_;
    std::vector<typename Voxel<M>::Ptr> voxels;
    const Parent& parent_;
};

template<typename M>
class DefaultVoxelDistFunctor
{
public:
    float dist(const typename Voxel<M>::Ptr& v1,const typename Voxel<M>::Ptr& v2,float seed_res)
    {
        float spatial_dist = arma::norm( v1->xyz() - v2->xyz() ) / seed_res;
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
    typedef std::function<float(const typename Voxel<M>::Ptr&,const typename Voxel<M>::Ptr&)> DistFunc;
    typedef __gnu_cxx::hash_multimap<uint64_t,arma::uvec> SuperVoxelsMap;
    typedef __gnu_cxx::hash_multimap<uint64_t,uint64_t> SuperVoxelAdjacency;
    typedef std::shared_ptr<M> MeshPtr;
    typedef unibn::Octree<arma::fvec,MeshOctreeContainer<M>> SuperVoxelOctree;
    typedef typename std::vector<typename Voxel<M>::Ptr>::iterator VoxelIter;
    typedef typename std::vector<typename SuperVoxel<M>::Ptr>::iterator SuperVoxelIter;
    SuperVoxelClustering(float v_res,float seed_res);
    SuperVoxelClustering(float v_res,float seed_res,bool);
    virtual ~SuperVoxelClustering();

    void setDistFunctor(DistFunc&);
    void input(MeshPtr mesh);
    void extract(arma::uvec&labels);
    void extract(SuperVoxelsMap&supervoxel_clusters);
    void getSupervoxelAdjacency(SuperVoxelAdjacency&label_adjacency)const;
protected:
    bool initCompute();
    void deinitCompute();
    void selectInitialSupervoxelSeeds(arma::uvec&indices);
    void createSupervoxels(arma::uvec& seed_indices);
    void expandSupervoxels(int depth);
    void makeSupervoxels(SuperVoxelsMap&);
protected:
    std::shared_ptr<SuperVoxelOctree> adjacency_octree_;
    std::vector<typename Voxel<M>::Ptr> voxels_;
    std::vector<typename SuperVoxel<M>::Ptr> supervoxels_;
    MeshPtr voxel_centroids_;
    MeshPtr input_;
private:
    float seed_resolution_;
    float voxel_resolution_;
    DistFunc dist_;
    friend class SuperVoxel<M>;
};
}
#include "supervoxelclustering.hpp"
#endif // SUPERVOXELCLUSTERING_H
