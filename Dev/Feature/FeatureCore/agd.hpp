#ifndef AGD_HPP
#define AGD_HPP
#include "agd.h"
#include <queue>
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
}
#endif // AGD_HPP
