#ifndef ROBUSTREGION_H
#define ROBUSTREGION_H
#include "normalizedcuts.h"
namespace Segmentation{
template<typename Mesh>
class RobustRegionDetection:public NormalizedCuts<Mesh>
{
public:
    using NormalizedCuts<Mesh>::configure;
    using NormalizedCuts<Mesh>::getLabel;
    using NormalizedCuts<Mesh>::computeW_Graph;
    using NormalizedCuts<Mesh>::decomposeGSP;
    using NormalizedCuts<Mesh>::clustering_Kmean;
    using NormalizedCuts<Mesh>::Y_;
    using NormalizedCuts<Mesh>::W_;
    using NormalizedCuts<Mesh>::gmm_;
public:
    RobustRegionDetection();
    void cutGraph(typename MeshBundle<Mesh>::Ptr m,arma::uvec&);
    arma::umat& base_segments(){return base_segments_;}
protected:
    void generate_base_segments();
    void solve_consensus_segment();
private:
    arma::uword n_base_segments_;
    arma::umat base_segments_;
};
}
#endif // ROBUSTREGION_H
