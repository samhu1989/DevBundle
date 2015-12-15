#include "sac_perpendicular_plane.h"
namespace Segmentation {
bool SAC_Perpendicular_Plane::isModelValid (const arma::vec& coeff)
{
    // Needs a valid model coefficients
    if (coeff.size () != 4)
    {
      std::cerr<<"Invalid number of model coefficients given "<<coeff.size ()<<std::endl;
      return (false);
    }

    // Check against template, if given
    if (eps_angle_ > 0.0)
    {
      double rad = arma::dot(coeff.head(3),axis_) / ( arma::norm(coeff.head(3))*arma::norm(axis_) );
      if (rad < -1.0) rad = -1.0;
      if (rad >  1.0) rad = 1.0;
      double angle_diff = std::abs(std::acos(rad));
      angle_diff = std::min(angle_diff, M_PI - angle_diff);
      // Check whether the current plane model satisfies our angle threshold criterion with respect to the given axis
      if (angle_diff > eps_angle_)
        return (false);
    }

    return (true);
}
uint64_t SAC_Perpendicular_Plane::countWithinDistance(arma::vec& coeff,double threshold)
{
    // Check if the model is valid given the user constraints
    if (!isModelValid (coeff))
    {
        std::cerr<<"Invalid coefficient"<<std::endl;
        return (0);
    }
    return SAC_Plane::countWithinDistance(coeff, threshold);
}
void SAC_Perpendicular_Plane::selectWithinDistance(arma::vec& coeff,double threshold,arma::uvec& inliers)
{    // Check if the model is valid given the user constraints
    if (!isModelValid (coeff))
    {
      std::cerr<<"Invalid coefficient"<<std::endl;
      inliers.reset();
      return;
    }
    SAC_Plane::selectWithinDistance(coeff,threshold,inliers);
}
}

