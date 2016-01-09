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
    int max_depth = 1 + static_cast<int> (1.8f*seed_resolution_/voxel_resolution_);
    expandSupervoxels (max_depth);
    makeLabels(labels);
    arma::uvec nz = arma::find(labels!=0);
    std::cerr<<"Numbers:"<<std::endl;
    std::cerr<<nz.size()<<std::endl;
    std::cerr<<input_->n_vertices()<<std::endl;
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
void SuperVoxelClustering<M>::getSupervoxelAdjacency(arma::Mat<uint16_t>&neighbors)
{
    SuperVoxelIter iter;
    std::vector<uint16_t> result;
    result.reserve(9*supervoxels_.size());
    arma::sp_umat mem(supervoxels_.size(),supervoxels_.size());
    __gnu_cxx::hash_map<uint32_t,uint32_t> labeltoindex;
    uint32_t index = 0;
    for(iter=supervoxels_.begin();iter!=supervoxels_.end();++iter)
    {
        SuperVoxel<M>& supervoxel = **iter;
        std::vector<uint32_t> neighbors;
        supervoxel.getNeighbors(neighbors);
        uint32_t label = supervoxel.getLabel();
        labeltoindex[label] = index;
        std::vector<uint32_t>::iterator labeliter;
        for(labeliter=neighbors.begin();labeliter!=neighbors.end();++labeliter)
        {
            uint32_t s,b;
            if(label<*labeliter){
                b = label;
                s = *labeliter;
            }else if(label==*labeliter)
            {
                continue;
            }else{
                b = *labeliter;
                s = label;
            }
            if(0==mem(s,b))
            {
                mem(s,b)=1;
                result.push_back(uint16_t(s));
                result.push_back(uint16_t(b));
            }
        }
        ++index;
    }
    std::vector<uint16_t>::iterator riter;
    for(riter=result.begin();riter!=result.end();++riter)
    {
        *riter = uint16_t(labeltoindex[uint32_t(*riter)]);
    }
    neighbors = arma::Mat<uint16_t>(result.data(),2,result.size()/2,true,true);
}

template<typename M>
void SuperVoxelClustering<M>::getCentroids(M&cmesh)
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
void SuperVoxelClustering<M>::getCentroids(arma::fmat& centers)
{
    SuperVoxelIter iter;
    centers = arma::fmat(3,supervoxels_.size(),arma::fill::ones);
    int idx = 0;
    for(iter=supervoxels_.begin();iter!=supervoxels_.end();++iter)
    {
        SuperVoxel<M>& supervoxel = **iter;
        Voxel<M>& center = *supervoxel.centeroid();
        centers.col(idx) = center.xyz();
        ++idx;
    }
}

template<typename M>
void SuperVoxelClustering<M>::getCentroidNormals(arma::fmat& normals)
{
    SuperVoxelIter iter;
    normals = arma::fmat(3,supervoxels_.size(),arma::fill::ones);
    int idx = 0;
    for(iter=supervoxels_.begin();iter!=supervoxels_.end();++iter)
    {
        SuperVoxel<M>& supervoxel = **iter;
        Voxel<M>& center = *supervoxel.centeroid();
        normals.col(idx) = center.normal();
        ++idx;
    }
}

template<typename M>
void SuperVoxelClustering<M>::getCentroidColors(arma::Mat<uint8_t>& colors)
{
    SuperVoxelIter iter;
    colors = arma::Mat<uint8_t>(3,supervoxels_.size(),arma::fill::ones);
    int idx = 0;
    for(iter=supervoxels_.begin();iter!=supervoxels_.end();++iter)
    {
        SuperVoxel<M>& supervoxel = **iter;
        Voxel<M>& center = *supervoxel.centeroid();
        colors.col(idx) = arma::conv_to<arma::Mat<uint8_t>>::from(center.color());
        ++idx;
    }
}

template<typename M>
void SuperVoxelClustering<M>::getSizes(arma::uvec& sizes)
{
    SuperVoxelIter iter;
    sizes = arma::uvec(supervoxels_.size(),arma::fill::zeros);
    int idx = 0;
    for(iter=supervoxels_.begin();iter!=supervoxels_.end();++iter)
    {
        SuperVoxel<M>& supervoxel = **iter;
        arma::uvec pointindices;
        supervoxel.getPointIndices(pointindices);
        sizes(idx) = pointindices.size();
        ++idx;
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
    return true;
}

template<typename M>
void SuperVoxelClustering<M>::computeVoxelData()
{
    voxels_.clear();
    voxels_.reserve(adjacency_octree_->getLeafNum());
    for(int i=0;i<adjacency_octree_->getLeafNum();++i)
    {
        arma::uvec indices;
        adjacency_octree_->getPointofLeafAt(i,indices);
        if(!indices.is_empty())voxels_.emplace_back(new Vox(*input_,indices));
        else{ std::cerr<<"encounter an empty voxel in octree"<<std::endl; }
    }
    std::cerr<<"Compute Voxel Data N:"<<voxels_.size()<<std::endl;
    OctreeVoxelAdjacency<SuperVoxelOctree> adjacency(*adjacency_octree_);
    adjacency.template computeNeighbor<Vox>( voxels_ );
//    int min_neighbor_n = 1000;
//    int max_neighbor_n = 0;
//    VoxelIter iter;
//    size_t num = 0;
//    for(iter=voxels_.begin();iter!=voxels_.end();++iter)
//    {
//        int n = (*iter)->neighbors_.size();
//        if(n<min_neighbor_n)min_neighbor_n=n;
//        if(n>max_neighbor_n)max_neighbor_n=n;
//        num += (*iter)->indices().size();
//    }
//    std::cerr<<"num points:"<<num<<std::endl;
//    std::cerr<<"Neighbor N("<<min_neighbor_n<<","<<max_neighbor_n<<")"<<std::endl;
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
    supervoxels_.clear();
    std::vector<uint32_t>::iterator iter;
    uint32_t label = 0;
    for(iter=seed_indices.begin();iter!=seed_indices.end();++iter)
    {
        supervoxels_.emplace_back(new SuperVoxel<M>(*this,label));
        voxels_[*iter]->dist2parent_ = 0.0;
        supervoxels_.back()->addVoxel(voxels_[*iter]);
        voxels_[*iter]->parent_ = supervoxels_.back();
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
        std::sort(supervoxels_.begin(), supervoxels_.end(),SuperVoxel<M>(*this,0));
        while(0==supervoxels_.back()->size())
        {
            supervoxels_.pop_back();
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
    labels.fill(0);
    arma::uword sindex = 0;
    for(iter=supervoxels_.begin();iter!=supervoxels_.end();++iter)
    {
        arma::uvec indices;
        (*iter)->getPointIndices(indices);
        labels.elem(indices).fill(sindex+1);
        ++sindex;
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
    centroid_->normal() = arma::normalise(centroid_->normal());
    for(iter=voxels.begin();iter!=voxels.end();++iter)
    {
        (*iter)->dist2parent_ = parent_.dist_(centroid_,(*iter));
    }
}

template<typename M>
void SuperVoxel<M>::expand()
{
    VoxelVec new_owned;
    new_owned.reserve(9*voxels.size());
    //For each leaf belonging to this supervoxel
    VoxelIter voxeliter;
    VoxelIter neighborIter;
    VoxelIter vbegin = voxels.begin ();
    for ( voxeliter = voxels.begin (); voxeliter != voxels.end (); ++voxeliter)
    {
        VoxelIter nbegin = (*voxeliter)->neighbors_.begin();
        VoxelIter nend = (*voxeliter)->neighbors_.end();

        for(neighborIter=nbegin;neighborIter!=nend;++neighborIter)
        { 
//            std::cerr<<"address:"<<(*neighborIter).get()<<std::endl;
//            std::cerr<<"count:"<<(*neighborIter).use_count()<<std::endl;
            //TODO this is a shortcut, really we should always recompute distance
            if(0!=(*neighborIter)->parent_.use_count()&&(*neighborIter)->parent_)
            {
                if((*neighborIter)->parent_.get()==this)continue;
            }
            //Compute distance to the neighbor
            float dist = parent_.dist_(centroid_,(*neighborIter));
            //If distance is less than previous, we remove it from its owner's list
            //and change the owner to this and distance (we *steal* it!)
            if (dist < (*neighborIter)->dist2parent_)
            {
                (*neighborIter)->dist2parent_ = dist;
                if (0!=(*neighborIter)->parent_.use_count()&&(*neighborIter)->parent_)
                {
                    (*neighborIter)->parent_ -> removeVoxel(*neighborIter);
                }
                (*neighborIter)->parent_ = voxels.front()->parent_;
                new_owned.push_back(*neighborIter);
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
    if(voxels.empty()){
        std::cerr<<"try to remove from an empty supervoxel"<<std::endl;
        return;
    }
    bool find = false;
    VoxelIter iter;
    if(ptr!=voxels.back())
    {
        for(iter=voxels.begin();iter!=voxels.end();++iter)
        {
            if(ptr==*iter)
            {
                (*iter)->parent_.reset();
                *iter = voxels.back();
                find = true;
                break;
            }
        }
    }else find = true;
    if(find){
        voxels.pop_back();
    }
    else {
        std::cerr<<"try to remove a voxel("<<ptr->Id_<<") with parent("<<ptr->parent_->getLabel()<<") but not recorded by supervoxel("<<label_<<")"<<std::endl;
        std::cerr<<"who has:"<<std::endl;
        for(iter=voxels.begin();iter!=voxels.end();++iter)
        {
            std::cerr<<(*iter)->Id_<<",";
        }
    }
}

template<typename M>
void SuperVoxel<M>::addVoxel(typename Voxel<M>::Ptr ptr)
{
    if(label_==40){
        std::cerr<<"40add"<<ptr->Id_<<std::endl;
    }
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
