#ifndef SEGMENTATIONBASE
#define SEGMENTATIONBASE
#include <armadillo>
namespace Segmentation
{
class SegmentationBase{
  public:
    virtual void extract(arma::uvec& labels)=0;
    virtual ~SegmentationBase(){}
};
}
#endif // SEGMENTATIONBASE

