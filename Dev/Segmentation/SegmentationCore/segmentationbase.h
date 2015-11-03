#ifndef SEGMENTATIONBASE
#define SEGMENTATIONBASE
namespace Segmentation
{
class SegmentationBase{
  public:
    virtual void getLabel(arma::uvec&)=0;
};
}
#endif // SEGMENTATIONBASE

