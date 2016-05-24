#include "extractpatchfeature.h"
#include "nanoflann.hpp"
using namespace nanoflann;
void extract_patch_expand(DefaultMesh&i,arma::uvec&indices,DefaultMesh&o,int k)
{
    MeshKDTreeInterface<DefaultMesh> pts(i);
    KDTreeSingleIndexAdaptor<
            L2_Simple_Adaptor<float,MeshKDTreeInterface<DefaultMesh>>,
            MeshKDTreeInterface<DefaultMesh>,
            3,arma::uword>
            kdtree(3,pts,KDTreeSingleIndexAdaptorParams(5));
    kdtree.buildIndex();
    arma::uvec expand_indices(k);
    arma::fvec dists(k);
    std::vector<arma::uword> expanded_vec = arma::conv_to<std::vector<arma::uword>>::from(indices);
    float* points = (float*)i.points();
//    std::cerr<<"searching:"<<std::endl;
    expanded_vec.reserve(5000);
    for(int idx=0;idx<indices.n_rows;++idx)
    {
        kdtree.knnSearch(&points[3*indices(idx)],k,expand_indices.memptr(),dists.memptr());
        for(int n=0;n<k;++n)expanded_vec.push_back(expand_indices(n));
    }
    arma::uvec expanded_indices(expanded_vec.data(),expanded_vec.size(),false,true);
    arma::uvec  unique_indices = arma::unique(expanded_indices);
    extractMesh<DefaultMesh,DefaultMesh>(i,unique_indices,o);
}

void extract_patch_expand(DefaultMesh&i,arma::uvec&indices,DefaultMesh&o,float r)
{
    MeshKDTreeInterface<DefaultMesh> pts(i);
    KDTreeSingleIndexAdaptor<
            L2_Simple_Adaptor<float,MeshKDTreeInterface<DefaultMesh>>,
            MeshKDTreeInterface<DefaultMesh>,
            3,arma::uword>
            kdtree(3,pts,KDTreeSingleIndexAdaptorParams(5));
    kdtree.buildIndex();
    std::vector<std::pair<arma::uword,float>> expand_indices;
    std::vector<std::pair<arma::uword,float>>::iterator iiter;
    std::vector<arma::uword> expanded_vec = arma::conv_to<std::vector<arma::uword>>::from(indices);
    float* points = (float*)i.points();
    arma::fmat pmat(points,3,i.n_vertices(),false,true);
    arma::fmat ppmat = pmat.cols(indices);
    arma::fvec min = arma::min(ppmat,1);
//    std::cerr<<"searching:"<<std::endl;
    expanded_vec.reserve(5000);
    for(int idx=0;idx<indices.n_rows;++idx)
    {
        kdtree.radiusSearch(&points[3*indices(idx)],r,expand_indices,SearchParams());
        for(iiter=expand_indices.begin();iiter!=expand_indices.end();++iiter)
        {
            expanded_vec.push_back(iiter->first);
        }
    }
    arma::uvec expanded_indices(expanded_vec.data(),expanded_vec.size(),false,true);
    arma::uvec  unique_indices = arma::unique(expanded_indices);
    extractMesh<DefaultMesh,DefaultMesh>(i,unique_indices,o);
}
