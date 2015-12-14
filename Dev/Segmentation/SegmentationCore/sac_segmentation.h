#ifndef SAC_SEGMENTATION_H
#define SAC_SEGMENTATION_H
#include "segmentationcore_global.h"
#include "SAC/sac_model.h"
#include "segmentationbase.h"
namespace Segmentation {
class SEGMENTATIONCORESHARED_EXPORT SegmentationRANSAC:public SegmentationBase
{
public:
    SegmentationRANSAC(SAC_Model::Model);
    inline void input(const arma::fmat& inputs){sac_model_->input(inputs);}
    inline void input(const arma::fmat& inputs,const arma::uvec& indices){
        sac_model_->input(inputs.cols(indices));
    }
    inline void extract(arma::uvec &labels);
    inline void get_model(arma::vec&);
protected:
    bool computeModel(void);
private:
    SAC_Model::Ptr sac_model_;
    arma::vec model_coefficients_;
    double threshold_;
    int max_iterations_;
    int iterations_;
    double probability_;
    arma::uvec inliers_;
    arma::uvec model_;
};
}
#endif // SAC_SEGMENTATION_H
