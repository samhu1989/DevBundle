#include "sac_segmentation.h"
#include "SAC/sac_plane.h"
#include "SAC/sac_parallel_plane.h"
#include "SAC/sac_perpendicular_plane.h"
namespace Segmentation{
SegmentationRANSAC::SegmentationRANSAC(SAC_Model::Model m):
    SegmentationBase(),
    probability_(0.99),
    iterations_(0),
    threshold_(0.05),
    max_iterations_ (1000)
{
    switch(m)
    {
    case SAC_Model::PLANE:
        sac_model_.reset(new SAC_Plane());
        break;
    case SAC_Model::PARALLEL_PLANE:
        sac_model_.reset(new SAC_Parallel_Plane());
        break;
    case SAC_Model::PERPENDICULLAR_PLANE:
        sac_model_.reset(new SAC_Perpendicular_Plane());
        break;
    }
}
bool SegmentationRANSAC::computeModel(void)
{
    if (threshold_ == std::numeric_limits<double>::max())
    {
      std::cerr<<"no threshold set"<<std::endl;
      return false;
    }

    iterations_ = 0;
    int n_best_inliers_count = -INT_MAX;
    double k = 1.0;

    arma::uvec selection;
    arma::vec model_coefficients;

    double log_probability  = std::log (1.0 - probability_);
    double one_over_indices = 1.0 / static_cast<double> (sac_model_->getInputSize());

    int n_inliers_count = 0;
    unsigned skipped_count = 0;
    // supress infinite loops by just allowing 10 x maximum allowed iterations for invalid model parameters!
    const unsigned max_skip = max_iterations_ * 10;
    std::cerr<<"Iterate"<<std::endl;
    while (iterations_ < k && skipped_count < max_skip)
    {
      std::cerr<<"Get X samples which satisfy the model criteria"<<std::endl;
      sac_model_->getSamples (iterations_, selection);

      if (selection.is_empty ())
      {
        std::cerr<<"No samples could be selected!"<<std::endl;
        break;
      }

      std::cerr<<"Search for inliers in the point cloud for the current plane model M"<<std::endl;
      if (!sac_model_->computeModel (selection, model_coefficients))
      {
        //++iterations_;
        ++skipped_count;
        continue;
      }

      // Select the inliers that are within threshold_ from the model
      //sac_model_->selectWithinDistance (model_coefficients, threshold_, inliers);
      //if (inliers.empty () && k > 1.0)
      //  continue;

      n_inliers_count = sac_model_->countWithinDistance (model_coefficients, threshold_);

      // Better match ?
      if (n_inliers_count > n_best_inliers_count)
      {
        n_best_inliers_count = n_inliers_count;

        std::cerr<<"Save the current model/inlier/coefficients selection as being the best so far"<<std::endl;
        model_              = selection;
        model_coefficients_ = model_coefficients;

        std::cerr<<"Compute the k parameter (k=log(z)/log(1-w^n))"<<std::endl;
        double w = static_cast<double> (n_best_inliers_count) * one_over_indices;
        double p_no_outliers = 1.0 - pow (w, static_cast<double> (selection.size ()));
        p_no_outliers = (std::max) (std::numeric_limits<double>::epsilon (), p_no_outliers);       // Avoid division by -Inf
        p_no_outliers = (std::min) (1.0 - std::numeric_limits<double>::epsilon (), p_no_outliers);   // Avoid division by 0.
        k = log_probability / log (p_no_outliers);
      }

      ++iterations_;
      if (iterations_ > max_iterations_)
      {
        break;
      }
    }

    if (model_.is_empty ())
    {
      inliers_.clear ();
      return (false);
    }
    std::cerr<<"Get the set of inliers that correspond to the best model found so far"<<std::endl;
    std::cerr<<"n_best_inliers_count:"<<n_best_inliers_count<<std::endl;
    sac_model_->selectWithinDistance (model_coefficients_, threshold_, inliers_);
    std::cerr<<"Inliers:"<<inliers_.size()<<std::endl;
    return (true);
}
}

