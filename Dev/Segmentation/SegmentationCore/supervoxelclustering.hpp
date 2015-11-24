#include "supervoxelclustering.h"
#include "common.h"
#include "nanoflann.hpp"
#include <utility>
namespace Segmentation{
template<typename M>
SuperVoxelClustering<M>::SuperVoxelClustering(float v_res,float seed_res)
    :voxel_resolution_(v_res),seed_resolution_(seed_res)
{

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
void SuperVoxelClustering<M>::input(M *mesh)
{
    input_ = mesh;
}

template<typename M>
void SuperVoxelClustering<M>::extract(arma::uvec&labels)
{
    bool segmentation_is_possible = initCompute ();
    if ( !segmentation_is_possible )
    {
        deinitCompute ();
        return;
    }
    std::vector<uint32_t> seed_indices;
    selectInitialSupervoxelSeeds (seed_indices);
    createSupervoxels(seed_indices);
    int max_depth = static_cast<int> (1.8f*seed_resolution_/voxel_resolution_);
    expandSupervoxels (max_depth);
    makeLabels(labels);
    deinitCompute ();
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
    std::vector<uint32_t> seed_indices;
    selectInitialSupervoxelSeeds (seed_indices);
    createSupervoxels(seed_indices);
    int max_depth = static_cast<int> (1.8f*seed_resolution_/voxel_resolution_);
    expandSupervoxels (max_depth);
    makeSupervoxels (supervoxel_clusters);
    deinitCompute ();
}

template<typename M>
void SuperVoxelClustering<M>::getSupervoxelAdjacency(SuperVoxelAdjacency&label_adjacency)
{
    SuperVoxelIter iter;
    for(iter=supervoxels_.begin();iter!=supervoxels_.end();++iter)
    {
        SuperVoxel<M>& supervoxel = **iter;
        std::vector<uint32_t> neighbors;
        supervoxel.getNeighbors(neighbors);
        uint32_t label = supervoxel.getLabel();
        std::vector<uint32_t>::iterator labeliter;
        for(labeliter=neighbors.begin();labeliter!=neighbors.end();++labeliter)
        {
            std::pair<uint32_t,uint32_t> p;
            p.first = label;
            p.second = *labeliter;
            label_adjacency.insert(p);
        }
    }
}

template<typename M>
void SuperVoxelClustering<M>::getCentroidMesh(M&cmesh)
{
    std::vector<typename M::VertexHandle> vhandle;
    SuperVoxelIter iter;
    for(iter=supervoxels_.begin();iter!=supervoxels_.end();++iter)
    {
        SuperVoxel<M>& supervoxel = **iter;
        Voxel<M>& center = *supervoxel.centeroid();
        arma::fvec& c = center.xyz();
        vhandle.push_back( cmesh.add_vertex( typename M::Point( c(0),c(1),c(2) ) ) );
    }
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
    computeVoxelData();
}

template<typename M>
void SuperVoxelClustering<M>::computeVoxelData()
{
    for(int i=0;i<adjacency_octree_->getLeafNum();++i)
    {
        arma::uvec indices;
        adjacency_octree_->getPointofLeafAt(i,indices);
        voxels_.emplace_back(new Vox(*input_,indices));
    }
    OctreeVoxelAdjacency<SuperVoxelOctree> adjacency(*adjacency_octree_);
    adjacency.template computeNeighbor<Vox>( voxels_ );
}

template<typename M>
void SuperVoxelClustering<M>::deinitCompute()
{
    ;
}

template<typename M>
void SuperVoxelClustering<M>::selectInitialSupervoxelSeeds(std::vector<uint32_t>&indices_vec)
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
}

template<typename M>
void SuperVoxelClustering<M>::createSupervoxels(std::vector<uint32_t> &seed_indices)
{
    std::vector<uint32_t>::iterator iter;
    uint32_t label = 0;
    for(iter=seed_indices.begin();iter!=seed_indices.end();++iter)
    {
        supervoxels_.push_back(std::make_shared<SuperVoxel<M>>(*this,label));
        supervoxels_.back()->addVoxel(voxels_[*iter]);
        supervoxels_.back()->updateCentroid();
        ++label;
    }
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
void SuperVoxelClustering<M>::makeSupervoxels(SuperVoxelsMap& supervoxelsmap)
{
    typename std::vector<typename SuperVoxel<M>::Ptr>::iterator iter;
    for(iter=supervoxels_.begin();iter!=supervoxels_.end();++iter)
    {
        std::vector<uint32_t> indices;
        (*iter)->getPointIndices(indices);
        std::vector<uint32_t>::iterator piter;
        for(piter=indices.begin();piter!=indices.end();++piter)
        {
            supervoxelsmap.insert(
                        std::make_pair(
                            (*iter)->getLabel(),
                            *piter
                            )
                        );
        }
    }
}

template<typename M>
void SuperVoxelClustering<M>::makeLabels(arma::uvec&labels)
{
    typename std::vector<typename SuperVoxel<M>::Ptr>::iterator iter;
    labels = arma::uvec(input_->n_vertices());
    labels.fill(std::numeric_limits<uint64_t>::max());
    for(iter=supervoxels_.begin();iter!=supervoxels_.end();++iter)
    {
        arma::uvec indices;
        (*iter)->getPointIndices(indices);
        arma::uword sindex = (*iter)->getLabel();
        labels.elem(indices).fill(sindex);
    }
}

template<typename M>
void SuperVoxel<M>::getPointIndices(std::vector<uint32_t>&indices)
{
    VoxelIter voxeliter;
    for ( voxeliter = voxels.begin (); voxeliter != voxels.end (); ++voxeliter)
    {
        const arma::uvec& vindices = (*voxeliter)->indices();
        for(int i=0;i<vindices.n_cols;++i)
        {
            indices.push_back(vindices(i));
        }
    }
}

template<typename M>
void SuperVoxel<M>::getPointIndices(arma::uvec&indices)
{
    VoxelIter voxeliter;
    for ( voxeliter = voxels.begin (); voxeliter != voxels.end (); ++voxeliter)
    {
        const arma::uvec& vindices = (*voxeliter)->indices();
        indices = arma::join_cols(vindices,indices);
    }
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
    centroid_->color() /= float(voxels.size());
    centroid_->xyz() /= float(voxels.size());
    centroid_->normal() /= float(voxels.size());
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

template<typename M>
void SuperVoxel<M>::getNeighbors(std::vector<uint32_t>&indices)
{
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
            if (neighbor_voxel.parent_){
                indices.push_back(neighbor_voxel.parent_->getLabel());
            }
        }
    }
}

}
