#ifndef REGIONGROWING_H
#define REGIONGROWING_H
#include "segmentationbase.h"
namespace Segmentation{
template<typename M>
class RegionGrowing:public SegmentationBase
{
public:
    RegionGrowing();
    void extract(arma::uvec&labels);
    virtual ~RegionGrowing();
    int getMinClusterSize(){return min_pts_per_cluster_;}
    void setMinClusterSize(int min_cluster_size){min_pts_per_cluster_=min_cluster_size;}
    int getMaxClusterSize(){return max_pts_per_cluster_;}
    void setMaxClusterSize(int max_cluster_size){max_pts_per_cluster_=max_cluster_size;}


protected:

private:
    int min_pts_per_cluster_;
    int max_pts_per_cluster_;

};

}
#endif // REGIONGROWING_H
