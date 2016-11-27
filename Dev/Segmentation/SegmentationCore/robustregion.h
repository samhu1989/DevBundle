#ifndef ROBUSTREGION_H
#define ROBUSTREGION_H
#include "normalizedcuts.h"
#include "rpeac.h"
namespace Segmentation{
template<typename Mesh>
class RobustRegionDetection:public NormalizedCuts<Mesh>
{
public:
    using NormalizedCuts<Mesh>::configure;
    using NormalizedCuts<Mesh>::getLabel;
    using NormalizedCuts<Mesh>::computeW_Graph;
    using NormalizedCuts<Mesh>::computeW_Image;
    using NormalizedCuts<Mesh>::decomposeGPS;
    using NormalizedCuts<Mesh>::clustering_Kmean;
    using NormalizedCuts<Mesh>::Y_;
    using NormalizedCuts<Mesh>::W_;
    using NormalizedCuts<Mesh>::gmm_;
public:
    RobustRegionDetection();
    bool configure(Config::Ptr);
    void cutGraph(typename MeshBundle<Mesh>::Ptr m,arma::uvec&);
    void generate_base_segments(typename MeshBundle<Mesh>::Ptr m,const arma::uvec& mask,bool re_use=false);
    void generate_base_segments(typename MeshBundle<Mesh>::Ptr m,bool re_use=false);
    void generate_base_segment(const QImage&,bool re_use=false);
    void solve_consensus_segment(typename MeshBundle<Mesh>::Ptr m,arma::uvec&);
    void solve_consensus_segment(const QImage&,arma::uvec&);
    arma::umat& base_segments(){return base_segments_;}
protected:
    void generate_base_segments(bool re_use_gps=false);
private:
    Clustering::RPEAC peac;
    arma::uword n_base_segments_;
    arma::umat base_segments_;
    arma::uvec masked_;
};
}
#endif // ROBUSTREGION_H
