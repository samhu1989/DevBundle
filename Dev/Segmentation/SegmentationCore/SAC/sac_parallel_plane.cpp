#include "sac_parallel_plane.h"
namespace Segmentation {
bool SAC_Parallel_Plane::isModelValid (const arma::vec& coeff)
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
      // Obtain the plane normal
      if (std::abs(arma::dot(coeff.head(3),arma::conv_to<arma::vec>::from(axis_))) > std::sin(eps_angle_) )
        return  (false);
    }
    return (true);
}
uint64_t SAC_Parallel_Plane::countWithinDistance(arma::vec& coeff,double threshold)
{
    // Needs a valid set of model coefficients
    if (coeff.size () != 4)
    {
      std::cerr<<"Invalid number of model coefficients given "<<coeff.size()<<std::endl;
      return (0);
    }

    int nr_p = 0;

    // Iterate through the 3d points and calculate the distances from them to the plane
    for (size_t i = 0; i < inputs_.n_cols; ++i)
    {
      // Calculate the distance from the point to the plane normal as the dot product
      // D = (P-A).N/|N|
      arma::vec pt(4);
      pt.head(3) = arma::conv_to<arma::vec>::from(inputs_.col(i));
      pt(3) = 1.0;
      if (std::abs (arma::dot(coeff,pt)) < threshold)
        nr_p++;
    }
    return (nr_p);
}
void SAC_Parallel_Plane::selectWithinDistance(arma::vec& coeff,double threshold,arma::uvec& inliers)
{
    // Check if the model is valid given the user constraints
    if (!isModelValid (coeff))
    {
      inliers.reset();
      return;
    }
    SAC_Plane::selectWithinDistance(coeff,threshold,inliers);
}
}

