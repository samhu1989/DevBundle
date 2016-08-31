#ifndef ROBUSTREGION_HPP
#define ROBUSTREGION_HPP
#include "robustregion.h"
namespace Segmentation{
template<typename Mesh>
RobustRegionDetection<Mesh>::RobustRegionDetection():NormalizedCuts<Mesh>()
{
    n_segments_ = 5;
}
template<typename Mesh>
bool RobustRegionDetection<Mesh>::configure(Config::Ptr)
{
    ;
}
template<typename Mesh>
void RobustRegionDetection<Mesh>::cutGraph(typename MeshBundle<Mesh>::Ptr m,arma::uvec& label)
{
    computeW_Graph(m);
    generate_base_segments();
    solve_consensus_segment();
    getLabel(label);
}
template<typename Mesh>
void RobustRegionDetection<Mesh>::generate_base_segments()
{
    assert(W_&&(0<W_.use_count()));
    base_segments_ = arma::umat(W_->n_rows,n_segments_,arma::fill::zeros);
    for(arma::uword i_segments_ = 0 ; i_segments_ < n_segments_ ; ++ i_segments_ )
    {
        decomposeGSP();
        clustering_Kmean();
        base_segments_.col(i_segments_) = (gmm_.assign(Y_.t(),arma::eucl_dist)).t();
    }
}
template<typename Mesh>
void RobustRegionDetection<Mesh>::solve_consensus_segment()
{
    ;
}
}
#endif // ROBUSTREGION_H
