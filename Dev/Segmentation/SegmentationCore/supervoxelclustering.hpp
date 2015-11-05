#include "supervoxelclustering.h"
#include "common.h"
#include "nanoflann.hpp"
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
    dist_ = vDist;
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
void SuperVoxelClustering<M>::extract(SuperVoxelsMap&supervoxel_clusters)
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
    adjacency_octree_.reset(new SuperVoxelOctree);
    if(!adjacency_octree_)return false;
    MeshOctreeContainer<M> input_container(*input_);
    adjacency_octree_->initialize(input_container,unibn::OctreeParams(float(voxel_resolution_)));
    voxel_centroids_.reset(new M);
    MeshOctreeContainer<M> center_container(*voxel_centroids_);
    adjacency_octree_->getLeafCenters(center_container);
}

template<typename M>
void SuperVoxelClustering<M>::deinitCompute()
{
    ;
}

template<typename M>
void SuperVoxelClustering<M>::selectInitialSupervoxelSeeds(arma::uvec&indices)
{
    SuperVoxelOctree seed_octree;
    MeshOctreeContainer<M> container(*voxel_centroids_);
    seed_octree.initialize(container,unibn::OctreeParams(float(seed_resolution_)));
    M seed_mesh;
    MeshOctreeContainer<M> seed_container(seed_mesh);
    int num_seeds = seed_octree.getLeafCenters(seed_container);

    MeshKDTreeInterface<M> points(*voxel_centroids_);
    nanoflann::KDTreeSingleIndexAdaptor<
            nanoflann::L2_Simple_Adaptor<float,MeshKDTreeInterface<M>>,
            MeshKDTreeInterface<M>,
            3,arma::uword>
            kdtree(3,points,nanoflann::KDTreeSingleIndexAdaptorParams(3));
    kdtree.buildIndex();
    if(num_seeds!=seed_mesh.n_vertices())std::logic_error("the seed_mesh size doesn't add up");
    arma::uvec seed_indices_orig(num_seeds,arma::fill::ones);
    arma::fvec seed_distance(num_seeds,arma::fill::zeros);
    float* seeds_ptr = (float*)seed_mesh.points();
    arma::uword* out_indices = (arma::uword*)seed_indices_orig.memptr();
    float* out_dist = (float*)seed_distance.memptr();
    for(int i = 0 ; i < seed_mesh.n_vertices() ; ++i )
    {
        kdtree.knnSearch(&(seeds_ptr[3*i]),1,&(out_indices[i]),&(out_dist[i]));
    }
    float search_radius = 0.5f*seed_resolution_;
    float min_points = 0.05f*(search_radius)*(search_radius)*M_PI/(voxel_resolution_*voxel_resolution_);
    std::vector<arma::uword> indices_vec;
    indices_vec.reserve(seed_indices_orig.size());
    float* points_ptr = (float*)voxel_centroids_->points();
    std::vector<std::pair<arma::uword,float>> neighbors;
    for(int i = 0 ; i < seed_indices_orig.size() ; ++i)
    {
        int num = kdtree.radiusSearch(
                    &(points_ptr[3*seed_indices_orig[i]]),
                    search_radius,
                    neighbors,nanoflann::SearchParams()
                    );
        if( num > min_points )
        {
           indices_vec.push_back(seed_indices_orig(i));
        }
    }
    indices = arma::uvec(indices_vec);
}

template<typename M>
void SuperVoxelClustering<M>::createSupervoxels(arma::uvec& seed_indices)
{
    ;
}

template<typename M>
void SuperVoxelClustering<M>::expandSupervoxels(int depth)
{
    SuperVoxelIter iter;
    for(int i = 0;i < depth ;++i)
    {
        for(iter=supervoxels_.begin();iter!=supervoxels_.end();++iter)
        {
            (*iter)->expand();
        }
        for(iter=supervoxels_.begin();iter!=supervoxels_.end();++iter)
        {
            if(0==(*iter)->size())
            {
                //erase empty supervoxel
                if((iter+1)==supervoxels_.end())supervoxels_.pop_back();
                else{
                    *iter = supervoxels_.back();
                    supervoxels_.pop_back();
                }
            }else{
                (*iter)->updateCentroid();
            }
        }
    }
}

template<typename M>
void SuperVoxelClustering<M>::makeSupervoxels(SuperVoxelsMap &)
{
    ;
}

template<typename M>
void SuperVoxel<M>::updateCentroid()
{
    centroid_.reset(new Voxel<M>(*parent_.input_));
    VoxelIter iter;
    for(iter=voxels.begin();iter!=voxels.end();++iter)
    {
        centroid_->color() += (*iter)->color();
        centroid_->normal() += (*iter)->normal();
        centroid_->xyz() += (*iter)->xyz();
    }
    centroid_->color() /= voxels.size();
    centroid_->xyz() /= voxels.size();
    centroid_->normal() /= voxels.size();
}

template<typename M>
void SuperVoxel<M>::expand()
{
    VoxelVec new_owned;
    new_owned.reserve(9*voxels.size());
    //For each leaf belonging to this supervoxel
    VoxelIter voxeliter;
    VoxelIter neighborIter;
    for ( voxeliter = voxels.begin (); voxeliter != voxels.end (); ++voxeliter)
    {
        //for each neighbor of the leaf
        for(neighborIter=(*voxeliter)->neighbors_.begin();neighborIter!=(*voxeliter)->neighbors_.end();++neighborIter)
        {
            //Get a reference to the data contained in the leaf
            Voxel<M>& neighbor_voxel = **neighborIter;
            //TODO this is a shortcut, really we should always recompute distance
            if(neighbor_voxel.parent_.get() == this)
                continue;
            //Compute distance to the neighbor
            float dist = parent_.dist_(centroid_,*neighborIter);
            //If distance is less than previous, we remove it from its owner's list
            //and change the owner to this and distance (we *steal* it!)
            if (dist < neighbor_voxel.dist2parent_)
            {
                neighbor_voxel.dist2parent_ = dist;
                if (neighbor_voxel.parent_.get() != this)
                {
                    if (neighbor_voxel.parent_)
                        (neighbor_voxel.parent_)->removeVoxel(*neighborIter);
                    neighbor_voxel.parent_.reset(this);
                    new_owned.push_back (*neighborIter);
                }
            }
        }
    }
    //Push all new owned onto the owned leaf set
    VoxelIter new_owned_itr;
    for ( new_owned_itr = new_owned.begin (); new_owned_itr!=new_owned.end (); ++new_owned_itr)
    {
        voxels.push_back(*new_owned_itr);
    }
}

template<typename M>
void SuperVoxel<M>::removeVoxel(typename Voxel<M>::Ptr ptr)
{
    if(ptr!=voxels.back())
    {
        VoxelIter iter;
        for(iter=voxels.begin();iter!=voxels.end();++iter)
        {
            if(ptr==*iter)
            {
                *iter = voxels.back();
                break;
            }
        }
    }
    voxels.pop_back();
}

template<typename M>
void SuperVoxel<M>::addVoxel(typename Voxel<M>::Ptr ptr)
{
    voxels.push_back(ptr);
}

}
