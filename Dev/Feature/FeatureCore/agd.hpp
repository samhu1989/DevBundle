#ifndef AGD_HPP
#define AGD_HPP
#include "agd.h"
#include <queue>
#include <algorithm>
#include <QMultiHash>
namespace Feature{
template<typename Mesh>
void AGD<Mesh>::extract(const Mesh&,arma::vec&)
{

}
template<typename Mesh>
void AGD<Mesh>::extract(const VoxelGraph<Mesh>& graph , arma::vec& agd)
{
    arma::uword N = graph.voxel_centers.n_cols;
    arma::mat gd(N,N,arma::fill::zeros);
    gd.fill(std::numeric_limits<float>::max());
    gd.diag().zeros();
    std::priority_queue<Pair> to_be_update;
    for(uint32_t index=0;index<graph.voxel_neighbors.n_cols;++index)
    {
        arma::uword i = graph.voxel_neighbors(0,index);
        arma::uword j = graph.voxel_neighbors(1,index);
        double dist = arma::norm(graph.voxel_centers.col(i) - graph.voxel_centers.col(j));
        to_be_update.emplace(i,j,dist);
        gd(i,j) = dist;
        gd(j,i) = dist;
    }
    while(!to_be_update.empty())
    {
        assert( to_be_update.size() < N*N*N );
        arma::uword i = to_be_update.top().first_;
        arma::uword j = to_be_update.top().second_;
        double dist = gd(i,j);
        to_be_update.pop();
        arma::uword update_index = 0;
        while(update_index<N)
        {
//            std::cerr<<"update_index:"<<update_index<<std::endl;
            if( (update_index==i) || (update_index==j) ){
                ++update_index;
                continue;
            }
            double gdi = gd(update_index,i);
            double gdj = gd(update_index,j);
            if( gdi < double(std::numeric_limits<float>::max()) )
            {
                if( ( gdi + dist ) < gdj ){
                    gd(update_index,j) = ( gdi + dist );
                    gd(j,update_index) = ( gdi + dist );
                    if( update_index < j )
                    {
                        to_be_update.emplace(update_index,j,gdi+dist);
                    }else{
                        to_be_update.emplace(j,update_index,gdi+dist);
                    }
                }
            }else if( gdj < double(std::numeric_limits<float>::max()) )
            {
                if( ( gdj + dist ) < gdi ){
                    gd(update_index,i) = ( gdj + dist );
                    gd(i,update_index) = ( gdj + dist );
                    if( update_index < i )
                    {
                        to_be_update.emplace(update_index,i,gdj+dist);
                    }else{
                        to_be_update.emplace(i,update_index,gdj+dist);
                    }
                }
            }
            ++update_index;
        }
    }
    if(gd.max() == std::numeric_limits<float>::max())
    {
        std::cerr<<"Unconnected parts exist"<<std::endl;
        gd.save("./debug/gd.mat.arma",arma::raw_ascii);
    }
    agd = arma::vec(gd.n_rows,arma::fill::zeros);
    #pragma omp parallel for
    for( uint32_t r = 0 ; r < gd.n_rows ; ++r )
    {
        arma::rowvec rr = gd.row(r);
        arma::uvec mi = arma::find(rr>0);
        if(!mi.empty())agd(r) = double(arma::mean(rr(mi)));
    }
}
template<typename Mesh>
void AGD<Mesh>::extract_slow(const VoxelGraph<Mesh>& graph , arma::vec& agd)
{
    arma::uword N = graph.voxel_centers.n_cols;
    arma::mat gd(N,N,arma::fill::zeros);
    gd.fill(std::numeric_limits<float>::max());
    gd.diag().zeros();
    std::queue<std::pair<uint16_t,uint16_t>> to_be_update;
    for(uint32_t index=0;index<graph.voxel_neighbors.n_cols;++index)
    {
        arma::uword i = graph.voxel_neighbors(0,index);
        arma::uword j = graph.voxel_neighbors(1,index);
        float dist = arma::norm(graph.voxel_centers.col(i) - graph.voxel_centers.col(j));
        gd(i,j) = dist;
        gd(j,i) = dist;
        to_be_update.push(std::make_pair<uint16_t,uint16_t>(i,j));
    }
    while(!to_be_update.empty())
    {
        arma::uword i = to_be_update.front().first;
        arma::uword j = to_be_update.front().second;
        float dist = gd(i,j);
        to_be_update.pop();
        for(arma::uword update_index=0 ; update_index < N; ++ update_index )
        {
            if( (update_index==i) || (update_index==j) )continue;
            float gdi = gd(update_index,i);
            float gdj = gd(update_index,j);
            if( gdi < std::numeric_limits<float>::max() )
            {
                if( ( gdi + dist ) < gdj ){
                    gd(update_index,j) = ( gdi + dist );
                    gd(j,update_index) = ( gdi + dist );
                    if( update_index < j )
                    {
                        to_be_update.push(std::make_pair<uint16_t,uint16_t>(update_index,j));
                    }else{
                        to_be_update.push(std::make_pair<uint16_t,uint16_t>(j,update_index));
                    }
                }
            }else if( gdj < std::numeric_limits<float>::max() )
            {
                if( ( gdj + dist ) < gdi ){
                    gd(update_index,i) = ( gdj + dist );
                    gd(i,update_index) = ( gdj + dist );
                    if( update_index < i )
                    {
                        to_be_update.push(std::make_pair<uint16_t,uint16_t>(update_index,i));
                    }else{
                        to_be_update.push(std::make_pair<uint16_t,uint16_t>(i,update_index));
                    }
                }
            }
        }
    }
//    std::cerr<<"gd:"<<std::endl;
//    std::cerr<<gd<<std::endl;
    if(gd.max() == std::numeric_limits<float>::max())
    {
        std::cerr<<"Unconnected parts exist"<<std::endl;
        gd.save("./debug/gd.mat.arma",arma::raw_ascii);
    }
    agd = arma::vec(gd.n_rows,arma::fill::zeros);
    #pragma omp parallel for
    for( uint32_t r = 0 ; r < gd.n_rows ; ++r )
    {
        arma::rowvec rr = gd.row(r);
        arma::uvec mi = arma::find(rr>0);
        if(!mi.empty())agd(r) = double(arma::mean(rr(mi)));
    }
}
template<typename Mesh>
void AGD<Mesh>::extract_hash(const VoxelGraph<Mesh>& graph,arma::vec& agd)
{
    arma::uword N = graph.voxel_centers.n_cols;
    arma::mat gd(N,N,arma::fill::zeros);
    gd.fill(std::numeric_limits<float>::max());
    gd.diag().zeros();
    QMultiHash<uint16_t,uint16_t> hash;
    for(uint32_t index=0;index<graph.voxel_neighbors.n_cols;++index)
    {
        arma::uword i = graph.voxel_neighbors(0,index);
        arma::uword j = graph.voxel_neighbors(1,index);
        float dist = arma::norm(graph.voxel_centers.col(i) - graph.voxel_centers.col(j));
        gd(i,j) = dist;
        gd(j,i) = dist;
        hash.insert(i,j);
        hash.insert(j,i);
    }
    std::queue<Pair> to_be_update;
    uint16_t i = hash.keys().first();
    uint16_t j = hash.values(i).first();
    if(i<j)to_be_update.emplace(i,j,gd(i,j));
    else if(j<i)to_be_update.emplace(j,i,gd(j,i));
    while(!to_be_update.empty())
    {
        arma::uword i = to_be_update.front().first_;
        arma::uword j = to_be_update.front().second_;
        double dist = gd(i,j);
        to_be_update.pop();
        //expand to i's neighbor
        for(QMultiHash<uint16_t,uint16_t>::iterator iter = hash.find(i);iter!=hash.end() && iter.key()==i ;++iter)
        {
            uint16_t& m = *iter;
            if(m==j)continue;
            double sumd = dist + gd(i,m);
            if( gd(j,m) == std::numeric_limits<float>::max() )
            {
                gd(j,m) = sumd;
                gd(m,j) = sumd;
                hash.insert(j,m);
                hash.insert(m,j);
                if(j<m)to_be_update.emplace(j,m,gd(j,m));
                else if(m<j)to_be_update.emplace(m,j,gd(m,j));
            }else if( gd(j,m) > sumd )
            {
                gd(j,m) = sumd;
                gd(m,j) = sumd;
                if(j<m)to_be_update.emplace(j,m,gd(j,m));
                else if(m<j)to_be_update.emplace(m,j,gd(m,j));
            }
        }
        //expand to j's neighbor
        for(QMultiHash<uint16_t,uint16_t>::iterator iter = hash.find(j);iter!=hash.end() && iter.key()==j;++iter)
        {
            uint16_t& m = *iter;
            if(m==i)continue;
            double sumd = dist + gd(j,m);
            if( gd(i,m) == std::numeric_limits<float>::max() )
            {
                gd(i,m) = sumd;
                gd(m,i) = sumd;
                hash.insert(i,m);
                hash.insert(m,i);
                if(i<m)to_be_update.emplace(i,m,gd(i,m));
                else if(m<i)to_be_update.emplace(m,i,gd(m,i));
            }else if( gd(i,m) > sumd )
            {
                gd(i,m) = sumd;
                gd(m,i) = sumd;
                if(i<m)to_be_update.emplace(i,m,gd(i,m));
                else if(m<i)to_be_update.emplace(m,i,gd(m,i));
            }
        }
    }
    if(gd.max() == std::numeric_limits<float>::max())
    {
        std::cerr<<"Unconnected parts exist"<<std::endl;
        gd.save("./debug/gd.mat.arma",arma::raw_ascii);
    }
    agd = arma::vec(gd.n_rows,arma::fill::zeros);
    #pragma omp parallel for
    for( uint32_t r = 0 ; r < gd.n_rows ; ++r )
    {
        arma::rowvec rr = gd.row(r);
        arma::uvec mi = arma::find(rr>0);
        if(!mi.empty())agd(r) = double(arma::mean(rr(mi)));
    }
}
}
#endif // AGD_HPP
