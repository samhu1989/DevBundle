#ifndef SEGMENTATIONBASE
#define SEGMENTATIONBASE
#include <armadillo>
namespace Segmentation
{
class SegmentationBase{
public:
    virtual void setIndices(arma::uvec& indices){indices_=indices;}
    virtual void extract(arma::uvec& labels)=0;
    virtual ~SegmentationBase(){}
protected:
    arma::uvec indices_;
};
}
#endif // SEGMENTATIONBASE

