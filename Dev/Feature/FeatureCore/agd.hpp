#ifndef AGD_HPP
#define AGD_HPP
#include "agd.h"
namespace Feature{
template<typename Mesh>
void AGD<Mesh>::extract(const Mesh&,arma::vec&)
{

}
template<typename Mesh>
void AGD<Mesh>::extract(const VoxelGraph<Mesh>& graph , arma::vec& agd)
{
    arma::uword N = graph.voxel_centers.n_cols;
    arma::fmat gd(N,N,arma::fill::zeros);
    gd.fill(std::numeric_limits<float>::max());
    gd.diag().zeros();
    for(uint32_t index=0;index<graph.voxel_neighbors.n_cols;++index)
    {
        arma::uword i = graph.voxel_neighbors(0,index);
        arma::uword j = graph.voxel_neighbors(1,index);
        float dist = arma::norm(graph.voxel_centers.col(i) - graph.voxel_centers.col(j));
        gd(i,j) = dist;
        gd(j,i) = dist;
        for(uint32_t update_index=0;update_index<N;++update_index)
        {
            if( (update_index==i) || (update_index==j) )continue;
            float gdi = gd(update_index,i);
            float gdj = gd(update_index,j);
            if( gdi < std::numeric_limits<float>::max() )
            {
                if( ( gdi + dist ) < gdj ){
                    gd(update_index,j) = ( gdi + dist );
                    gd(j,update_index) = ( gdi + dist );
                }
            }else if( gdj < std::numeric_limits<float>::max() )
            {
                if( ( gdj + dist ) < gdi ){
                    gd(update_index,i) = ( gdj + dist );
                    gd(i,update_index) = ( gdj + dist );
                }
            }
        }
    }
    agd = arma::conv_to<arma::vec>::from(arma::mean(gd,1));
}
}
#endif // AGD_HPP
