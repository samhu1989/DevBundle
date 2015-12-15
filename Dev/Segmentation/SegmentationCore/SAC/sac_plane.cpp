#include "sac_plane.h"
namespace Segmentation {
bool SAC_Plane::isSampleGood (const std::vector<int> &samples)const
{
    // Need an extra check in case the sample selection is empty
    if (samples.empty ())
      return (false);
    // Get the values at the two points
    arma::fvec p0 = inputs_.col(samples[0]);
    arma::fvec p1 = inputs_.col(samples[1]);
    arma::fvec p2 = inputs_.col(samples[2]);
    arma::fvec dy1dy2 = (p1-p0) / (p2-p0);
    return ( (dy1dy2[0] != dy1dy2[1]) || (dy1dy2[2] != dy1dy2[1]) );
}
bool SAC_Plane::isModelValid (const arma::vec& coeff)
{
    if (coeff.size() != 4)
    {
      std::cerr<<"Invalid number of model coefficients given "<<coeff.size()<<std::endl;
      return (false);
    }
    return (true);
}
bool SAC_Plane::computeModel(arma::uvec&inliers,arma::vec& coeff)
{
    if( inliers.size() != 3)
    {
        std::cerr<<"Invalid set of samples given"<<inliers.size()<<std::endl;
        return (false);
    }

    arma::fvec p0 = inputs_.col(inliers(0));
    arma::fvec p1 = inputs_.col(inliers(1));
    arma::fvec p2 = inputs_.col(inliers(2));

    // Compute the segment values (in 3d) between p1 and p0
    arma::vec p1p0 = arma::conv_to<arma::vec>::from(p1 - p0);
    // Compute the segment values (in 3d) between p2 and p0
    arma::vec p2p0 = arma::conv_to<arma::vec>::from(p2 - p0);

    // Avoid some crashes by checking for collinearity here
    arma::vec dy1dy2 = p1p0 / p2p0;
    if ( (dy1dy2[0] == dy1dy2[1]) && (dy1dy2[2] == dy1dy2[1]) )          // Check for collinearity
        return (false);

    // Compute the plane coefficients from the 3 given points in a straightforward manner
    // calculate the plane normal n = (p2-p1) x (p3-p1) = cross (p2-p1, p3-p1)
    coeff.resize (4);
    coeff[0] = p1p0[1] * p2p0[2] - p1p0[2] * p2p0[1];
    coeff[1] = p1p0[2] * p2p0[0] - p1p0[0] * p2p0[2];
    coeff[2] = p1p0[0] * p2p0[1] - p1p0[1] * p2p0[0];
    coeff[3] = 0;

    // Normalize
    coeff.head(3) = arma::normalise(coeff.head(3));

    // ... + d = 0
    coeff[3] = -1 * (arma::dot(coeff.head(3),arma::conv_to<arma::vec>::from(p0)));

    if(!isModelValid(coeff))
    {
        return false;
    }
    return (true);
}
void SAC_Plane::optimizeModel(const arma::uvec&inliers,
                           const arma::vec& coeff,
                           arma::vec& optimized_coeff)
{
    std::cerr<<"Needs a valid set of model coefficients"<<std::endl;
    if (coeff.size () != 4)
    {
      std::cerr<<"Invalid number of model coefficients given "<<coeff.size ()<<std::endl;
      optimized_coeff = coeff;
      return ;
    }

    std::cerr<<"Need at least 3 points to estimate a plane"<<std::endl;
    if (inliers.size() < 3)
    {
      std::cerr<<"Not enough inliers found to support a model"<<inliers.size ()<<"! Returning the same coefficients.\n"<<std::endl;
      optimized_coeff = coeff;
      return;
    }

    arma::fvec plane_normal;
    float plane_dist;
    float curvature;
    std::cerr<<"Use Least-Squares to fit the plane through all the given sample points and find out its coefficients"<<std::endl;

    arma::fmat neighborhood = inputs_.cols(inliers);
    arma::fvec xyz_centroid = arma::mean(neighborhood,1);
    fitPlane(xyz_centroid,neighborhood,plane_normal,curvature,plane_dist);

    std::cerr<<"Hessian form (D = nc . p_plane (centroid here) + p)"<<std::endl;
    optimized_coeff.resize (4);
    optimized_coeff[0] = plane_normal [0];
    optimized_coeff[1] = plane_normal [1];
    optimized_coeff[2] = plane_normal [2];
    optimized_coeff[3] = plane_dist;

    std::cerr<<"Make sure it results in a valid model"<<std::endl;
    if (!isModelValid (optimized_coeff))
    {
      optimized_coeff = coeff;
    }
}
uint64_t SAC_Plane::countWithinDistance(arma::vec& coeff,double threshold)
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
void SAC_Plane::selectWithinDistance(arma::vec& coeff,double threshold,arma::uvec& inliers)
{
    // Needs a valid set of model coefficients
    if (coeff.size () != 4)
    {
      std::cerr<<"Invalid number of model coefficients given "<<coeff.size ()<<std::endl;
      return;
    }

    int nr_p = 0;
    std::vector<arma::uword> inliersvec;
    inliersvec.reserve(inputs_.n_cols);
    error_sqr_dists_.clear();
    error_sqr_dists_.reserve(inputs_.n_cols);

    // Iterate through the 3d points and calculate the distances from them to the plane
    for (size_t i = 0; i < inputs_.n_cols; ++i)
    {
      // Calculate the distance from the point to the plane normal as the dot product
      // D = (P-A).N/|N|
        arma::vec pt(4);
        pt.head(3) = arma::conv_to<arma::vec>::from(inputs_.col(i));
        pt(3) = 1.0;

      float distance = std::abs (arma::dot(coeff,pt));

      if (distance < threshold)
      {
        // Returns the indices of the points whose distances are smaller than the threshold
        inliersvec.push_back(i);
        error_sqr_dists_.push_back(distance);
        ++nr_p;
      }
    }
    inliers = arma::conv_to<arma::uvec>::from(inliersvec);
}
SAC_Model::Model SAC_Plane::type()
{
    return SAC_Model::PLANE;
}
}

