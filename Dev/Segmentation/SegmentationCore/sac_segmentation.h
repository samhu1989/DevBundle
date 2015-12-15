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
    inline void input(const arma::fmat& inputs)
    {
        iterations_ = 0;
        sac_model_->input(inputs);
        input_full_size_ = inputs.n_cols;
        input_indices_.reset();
    }
    inline void input(const arma::fmat& inputs,const arma::uvec& indices)
    {
        iterations_ = 0;
        sac_model_->input(inputs.cols(indices));
        input_full_size_ = inputs.n_cols;
        input_indices_   = indices;
    }
    inline void extract(arma::uvec &labels)
    {
        computeModel();
        sac_model_->optimizeModel(inliers_,model_coefficients_,model_coefficients_);
        std::cerr<<model_coefficients_<<std::endl;
        sac_model_->selectWithinDistance(model_coefficients_,threshold_,inliers_);
        getInliersByLabel(labels);

    }
    inline void get_model(arma::vec& coeff){coeff=model_coefficients_;}
    inline void getInliersByLabel(arma::uvec& label)
    {
        if(label.is_empty())label = arma::uvec(input_full_size_,arma::fill::ones);
        if(input_indices_.is_empty())label(inliers_).fill(0);
        else label(input_indices_(inliers_)).fill(0);
    }
    inline void setThreshold(const float& th){threshold_=th;}
    inline void setAxis(const arma::fvec&axis){sac_model_->setAxis(axis);}
    inline void setEpsAngle(const float&eps){sac_model_->setEpsAngle(eps);}
protected:
    bool computeModel(void);
private:
    SAC_Model::Ptr sac_model_;
    uint64_t input_full_size_;
    arma::uvec input_indices_;
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
