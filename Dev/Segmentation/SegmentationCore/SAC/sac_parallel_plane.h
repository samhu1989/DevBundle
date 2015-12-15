#ifndef SAC_PARALLEL_PLANE_H
#define SAC_PARALLEL_PLANE_H
#include "sac_plane.h"
#include <memory>
namespace Segmentation {
class SAC_Parallel_Plane:public SAC_Plane
{
public:
    SAC_Parallel_Plane():SAC_Plane(),axis_(3,arma::fill::zeros),eps_angle_(M_PI/180.0){}
    virtual inline void setAxis(const arma::fvec&axis){axis_=axis;}
    virtual inline void setEpsAngle(const float&eps){eps_angle_=eps;}
    virtual uint64_t countWithinDistance(arma::vec& coeff,double threshold);
    virtual void selectWithinDistance(arma::vec& coeff,double threshold,arma::uvec& inliers);
    virtual Model type(){return PARALLEL_PLANE;}
protected:
    virtual inline bool isModelValid (const arma::vec& coeff);
private:
    std::vector<float> error_sqr_dists_;
    arma::fvec axis_;
    float eps_angle_;
};
}
#endif // SAC_PARALLEL_PLANE_H
