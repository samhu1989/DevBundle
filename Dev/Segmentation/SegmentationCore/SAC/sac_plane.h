#ifndef SAC_PLANE_H
#define SAC_PLANE_H
#include "sac_model.h"
#include "segmentationcore_global.h"
#include <memory>
namespace Segmentation {
class SEGMENTATIONCORESHARED_EXPORT SAC_Plane:public SAC_Model
{
public:
    using SAC_Model::inputs_;
    using SAC_Model::shuffled_indices_;
    SAC_Plane():SAC_Model(){}
    virtual bool computeModel(arma::uvec&inliers,arma::vec& coeff);
    virtual void optimizeModel(const arma::uvec&inliers,
                               const arma::vec& coeff,
                               arma::vec& optimized_coeff);
    virtual uint64_t countWithinDistance(arma::vec& coeff,double threshold);
    virtual void selectWithinDistance(arma::vec& coeff,double threshold,arma::uvec& inliers);
    virtual Model type();
protected:
    virtual bool isSampleGood(const std::vector<int> &samples)const;
    virtual inline bool isModelValid (const arma::vec& coeff);
    virtual int  getSampleSize()const{return 3;}
private:
    std::vector<float> error_sqr_dists_;
};
}
#endif
