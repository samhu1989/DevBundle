#include "supervoxelclustering.h"
#include "common.h"
namespace Segmentation{
template<typename M>
SuperVoxelClustering<M>::SuperVoxelClustering(float v_res,float seed_res)
{
    ;
}
template<typename M>
SuperVoxelClustering<M>::SuperVoxelClustering(float v_res,float seed_res,bool)
{
    ;
}

template<typename M>
SuperVoxelClustering<M>::~SuperVoxelClustering()
{
    ;
}

template<typename M>
void SuperVoxelClustering<M>::setDistFunctor(DistFunc& vDist)
{
    vDist_ = vDist;
}

template<typename M>
void SuperVoxelClustering<M>::input(MeshPtr mesh)
{
    ;
}

template<typename M>
void SuperVoxelClustering<M>::extract(arma::uvec&labels)
{
    ;
}

template<typename M>
void SuperVoxelClustering<M>::extract(SuperVoxels&supervoxel_clusters)
{
    bool segmentation_is_possible = initCompute ();
    if ( !segmentation_is_possible )
    {
        deinitCompute ();
        return;
    }

    arma::uvec seed_indices;
    selectInitialSupervoxelSeeds (seed_indices);

    createSupervoxels(seed_indices);

    int max_depth = static_cast<int> (1.8f*seed_resolution_/voxel_resolution_);
    expandSupervoxels (max_depth);

    makeSupervoxels (supervoxel_clusters);

    deinitCompute ();
}

template<typename M>
void SuperVoxelClustering<M>::getSupervoxelAdjacency(SuperVoxelAdjacency&label_adjacency)const
{
    ;
}

template<typename M>
bool SuperVoxelClustering<M>::initCompute()
{
    if(!input_)return false;
}

template<typename M>
void SuperVoxelClustering<M>::deinitCompute()
{
    ;
}

template<typename M>
void SuperVoxelClustering<M>::selectInitialSupervoxelSeeds(arma::uvec&indices)
{
    unibn::Octree<arma::fvec,MeshOctreeContainer<DefaultMesh>> seed_octree;
}

template<typename M>
void SuperVoxelClustering<M>::createSupervoxels(arma::uvec& seed_indices)
{

}

template<typename M>
void SuperVoxelClustering<M>::expandSupervoxels(int depth)
{
    ;
}

template<typename M>
void SuperVoxelClustering<M>::makeSupervoxels(SuperVoxels&)
{
    ;
}

}
