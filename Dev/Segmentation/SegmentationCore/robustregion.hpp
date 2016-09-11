#ifndef ROBUSTREGION_HPP
#define ROBUSTREGION_HPP
#include "robustregion.h"

namespace Segmentation{
template<typename Mesh>
RobustRegionDetection<Mesh>::RobustRegionDetection():NormalizedCuts<Mesh>()
{
    n_base_segments_ = 5;
}
template<typename Mesh>
void RobustRegionDetection<Mesh>::cutGraph(typename MeshBundle<Mesh>::Ptr m,arma::uvec& label)
{
    computeW_Graph(m);
    generate_base_segments();
    solve_consensus_segment(m,label);
}
template<typename Mesh>
bool RobustRegionDetection<Mesh>::configure(Config::Ptr config)
{
    if(!NormalizedCuts<Mesh>::configure(config))return false;
    if(config->has("Base_Seg_Num")){
        n_base_segments_ = config->getInt("Base_Seg_Num");
    }
    if(!peac.configure(config))return false;
    return true;
}
template<typename Mesh>
void RobustRegionDetection<Mesh>::generate_base_segments(typename MeshBundle<Mesh>::Ptr m)
{
    computeW_Graph(m);
    generate_base_segments();
}

template<typename Mesh>
void RobustRegionDetection<Mesh>::solve_consensus_segment(typename MeshBundle<Mesh>::Ptr m, arma::uvec& label)
{
    arma::umat edge = arma::conv_to<arma::umat>::from(m->graph_.voxel_neighbors);
    assert(base_segments_.n_rows==m->graph_.voxel_centers.n_cols);
    peac.compute(base_segments_.t(),edge,label);
}

template<typename Mesh>
void RobustRegionDetection<Mesh>::generate_base_segments()
{
    assert(W_&&(0<W_.use_count()));
    if( base_segments_.n_rows!=W_->n_rows || base_segments_.n_cols!=n_base_segments_)
    {
        base_segments_ = arma::umat(W_->n_rows,n_base_segments_,arma::fill::zeros);
    }
    for(arma::uword i_segments_ = 0 ; i_segments_ < n_base_segments_ ; ++ i_segments_ )
    {
        decomposeGPS();
        clustering_Kmean();
        base_segments_.col(i_segments_) = (gmm_.assign(Y_.t(),arma::eucl_dist)).t();
    }
}
}
#endif // ROBUSTREGION_H
