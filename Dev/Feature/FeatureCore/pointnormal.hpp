#include "pointnormal.h"
#include "nanoflann.hpp"
#include "featurecore.h"
using namespace nanoflann;
namespace Feature{
template<typename M>
void computePointNormal(M& mesh,float r,int k)
{
    MeshKDTreeInterface<M> points(mesh);
    KDTreeSingleIndexAdaptor<
            L2_Simple_Adaptor<float,MeshKDTreeInterface<M>>,
            MeshKDTreeInterface<M>,
            3,arma::uword>
            kdtree(3,points,KDTreeSingleIndexAdaptorParams(9));
    kdtree.buildIndex();
    if(!mesh.has_vertex_normals())
    {
        mesh.request_vertex_normals();
    }
    nanoflann::SearchParams param;
    float* point_ptr = (float*)mesh.points();
    float* normal_ptr = (float*)mesh.vertex_normals();
    arma::fmat X(point_ptr,3,mesh.n_vertices(),false,true);
    int idx;
    if( k < 4 )k = 4;
    for(idx=0;idx<mesh.n_vertices();++idx)
    {
        std::vector<std::pair<arma::uword,float>> neighbor;
        std::vector<std::pair<arma::uword,float>>::iterator iter;
        arma::fmat neighborMat;
        int neighborSize;
        if(r!=0.0){
            neighborSize = kdtree.radiusSearch(point_ptr,r,neighbor,param);
            arma::uvec indices(neighbor.size());
            int indicesidx = 0;
            for (iter==neighbor.begin();iter!=neighbor.end();++iter) {
                indices(indicesidx) = iter->first;
                ++indicesidx;
            }
            neighborMat = X.cols(indices);
        }
        else
        {
            arma::uword index[k];
            float dist[k];
            kdtree.knnSearch(point_ptr,k,&index[0],&dist[0]);
            arma::uvec indices(index,k,false,true);
            neighborMat = X.cols(indices);
        }
        arma::fvec center(point_ptr,3,false,true);
        arma::fvec normal(normal_ptr,3,false,true);
        if(neighborMat.n_cols<4)
        {
            normal.fill(std::numeric_limits<float>::quiet_NaN());
        }else{
            float d2o;//distance to origin
            fitPlane(center,neighborMat,normal,d2o);
        }
        point_ptr+=3;
        normal_ptr+=3;
    }
}
template<typename M>
void computePointNormal(M& mesh,std::shared_ptr<float>& curvature,float r,int k)
{
    MeshKDTreeInterface<M> points(mesh);
    KDTreeSingleIndexAdaptor<
            L2_Simple_Adaptor<float,MeshKDTreeInterface<M>>,
            MeshKDTreeInterface<M>,
            3,arma::uword>
            kdtree(3,points,KDTreeSingleIndexAdaptorParams(9));
    kdtree.buildIndex();
    if(!mesh.has_vertex_normals())
    {
        mesh.request_vertex_normals();
    }
    nanoflann::SearchParams param;
    float* point_ptr = (float*)mesh.points();
    float* normal_ptr = (float*)mesh.vertex_normals();
    if(!curvature)curvature.reset(new float[mesh.n_vertices()]);
    float* c_ptr = curvature.get();
    arma::fmat X(point_ptr,3,mesh.n_vertices(),false,true);
    int idx;
    if( k < 4 )k = 4;
    for(idx=0;idx<mesh.n_vertices();++idx)
    {
        std::vector<std::pair<arma::uword,float>> neighbor;
        std::vector<std::pair<arma::uword,float>>::iterator iter;
        arma::fmat neighborMat;
        int neighborSize;
        if(r!=0.0){
            neighborSize = kdtree.radiusSearch(point_ptr,r,neighbor,param);
            arma::uvec indices(neighbor.size());
            int indicesidx = 0;
            for (iter==neighbor.begin();iter!=neighbor.end();++iter) {
                indices(indicesidx) = iter->first;
                ++indicesidx;
            }
            neighborMat = X.cols(indices);
        }
        else
        {
            arma::uword index[k];
            float dist[k];
            kdtree.knnSearch(point_ptr,k,&index[0],&dist[0]);
            arma::uvec indices(index,k,false,true);
            neighborMat = X.cols(indices);
        }
        arma::fvec center(point_ptr,3,false,true);
        arma::fvec normal(normal_ptr,3,false,true);
        if(neighborMat.n_cols<4)
        {
            normal.fill(std::numeric_limits<float>::quiet_NaN());
        }else{
            float d2o;//distance to origin
            fitPlane(center,neighborMat,normal,c_ptr[idx],d2o);
            //flip normal towards (0,0,0)
            if(arma::dot(normal,center)>0)
            {
                normal*=-1.0;
            }
        }
        point_ptr+=3;
        normal_ptr+=3;
    }
}
}

