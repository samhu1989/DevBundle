#ifndef GDCOORD_HPP
#define GDCOORD_HPP
#include "gdcoord.h"
#include <armadillo>
#include <queue>
namespace Feature{
struct Node
{
    float* dist_;
    arma::uword id_;
    bool visited_;
    friend bool operator<(Node a, Node b)
    {
        return *a.dist_ > *b.dist_;
    }
};
template<typename Mesh>
void GDCoord<Mesh>::extract(const VoxelGraph<Mesh>& graph, const arma::uvec& axis_index ,arma::fmat& feature)
{
    std::cerr<<"GDCoord<Mesh>::extract()"<<std::endl;
    if(feature.n_rows!=axis_index.size()||feature.n_cols!=graph.size())feature = arma::fmat(axis_index.size(),graph.size());
    feature.fill( std::numeric_limits<float>::max() - 1.0 );
    arma::uvec axis = axis_index - 1;
    //build connected matrix
    size_t N = graph.size();
    arma::sp_fmat w(N,N);
    for(uint32_t index=0;index < graph.voxel_neighbors.n_cols;++index)
    {
        arma::uword i = graph.voxel_neighbors(0,index);
        arma::uword j = graph.voxel_neighbors(1,index);
        float dist = arma::norm(graph.voxel_centers.col(i) - graph.voxel_centers.col(j));
        w(i,j) = dist;
        w(j,i) = dist;
    }
    arma::uword dim = axis.size();
    std::cerr<<"dim:"<<dim<<std::endl;
    float* f_ptr_ = feature.memptr();
    #pragma omp parallel for
    for( int i = 0; i < dim ; ++i )
    {
        //Dijkstra
        std::vector<Node> d(N);
        std::priority_queue<Node> q;
        for(int j = 0; j < N; j++)
        {
            Node& node = d[j];
            node.id_ = j;
            node.dist_ = &f_ptr_[j*dim+i];
            node.visited_ = false;
        }
        q.push(d[axis(i)]);
        *d[axis(i)].dist_  = 0.0;
        while(!q.empty())
        {
            Node cd = q.top();
            q.pop();
            int u = cd.id_;
            if(cd.visited_)
                continue;
            cd.visited_ = true;
            for(int v = 0; v < N; v++)
            {
                if(v != u && !d[v].visited_ && w(u,v) != 0 && (*d[v].dist_) > ( *d[u].dist_ + w(u,v) ))
                {
                    (*d[v].dist_) = (*d[u].dist_ + w(u,v));
                    q.push(d[v]);
                }
            }
        }
    }
//    std::cerr<<"f:"<<feature<<std::endl;
    feature += 1.0 ;
    feature = 1.0 / feature;
    if(!feature.is_finite())
    {
        std::cerr<<"infinite feature GDCoord<Mesh>::extract"<<std::endl;
    }
}
}
#endif // GDCOORD_HPP
