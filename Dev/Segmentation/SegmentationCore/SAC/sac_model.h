#ifndef SAC_MODEL_H
#define SAC_MODEL_H
#include "segmentationcore_global.h"
#include <memory>
#include <armadillo>
#include <random>
#include <limits>
#include "nanoflann.hpp"
#include "common.h"
namespace Segmentation {
using namespace nanoflann;
class SEGMENTATIONCORESHARED_EXPORT SAC_Model
{
public:
    typedef enum{
        PLANE,
        PARALLEL_PLANE,
        PERPENDICULLAR_PLANE
    }Model;
    typedef std::shared_ptr<SAC_Model> Ptr;
    typedef KDTreeSingleIndexAdaptor<
            L2_Simple_Adaptor<float,ArmaKDTreeInterface<arma::fmat>>,
            ArmaKDTreeInterface<arma::fmat>,
            3,arma::uword
            > Search;
    typedef std::shared_ptr<Search> SearchPtr;
    SAC_Model():
        rnd_gen_(rnd_d_()),
        rnd_dist_(0,std::numeric_limits<int>::max()),
        samples_radius_(0.0)
    {
    }
    virtual ~SAC_Model(){}
    inline void input(const arma::fmat& inputs)
    {
        inputs_ = inputs;
        arma::uvec indices = arma::linspace<arma::uvec>(0,inputs_.n_cols-1,inputs_.n_cols);
        shuffled_indices_ = arma::shuffle(indices);
    }
    virtual inline void setAxis(const arma::fvec&){std::cerr<<"not implemented setAxis"<<std::endl;}
    virtual inline void setEpsAngle(const float&){std::cerr<<"not implemented setEpsAngle(float&)"<<std::endl;}
    inline uint64_t getInputSize(){
        return inputs_.n_cols;
    }
    inline void
    drawIndexSample (std::vector<int> &sample)
    {
//      std::cerr<<"drawIndexSample(std::vector<int> &sample)"<<std::endl;
      size_t sample_size = sample.size ();
      size_t index_size = shuffled_indices_.size ();
//      std::cerr<<"sample_size:"<<sample_size<<std::endl;
//      std::cerr<<"index_size:"<<index_size<<std::endl;
      for (unsigned int i = 0; i < sample_size; ++i)
        shuffled_indices_.swap_rows(i,i + ( rnd() % (index_size - i)));
//      std::cerr<<"copy to sample"<<std::endl;
        std::copy (shuffled_indices_.begin(), shuffled_indices_.begin()+sample_size, sample.begin() );
    }
    inline void
    drawIndexSampleRadius (std::vector<int> &sample)
    {
      size_t sample_size = sample.size ();
      size_t index_size = shuffled_indices_.size ();

      std::swap (shuffled_indices_[0], shuffled_indices_[0 + (rnd () % (index_size - 0))]);
      //const PointT& pt0 = (*input_)[shuffled_indices_[0]];

      std::vector<std::pair<arma::uword,float>> results;
      float* pts = (float*)inputs_.memptr();
      // If indices have been set when the search object was constructed,
      // radiusSearch() expects an index into the indices vector as its
      // first parameter. This can't be determined efficiently, so we use
      // the point instead of the index.
      // Returned indices are converted automatically.
      if(samples_radius_search_&&0!=samples_radius_search_.use_count())
          samples_radius_search_->radiusSearch(&pts[shuffled_indices_[0]],
                                            samples_radius_, results,SearchParams());
      else{
          drawIndexSample(sample);
          return;
      }

      if ( results.size() < sample_size - 1)
      {
        // radius search failed, make an invalid model
        for(unsigned int i = 1; i < sample_size; ++i)
          shuffled_indices_[i] = shuffled_indices_[0];
      }
      else
      {
        arma::uvec indices = arma::linspace<arma::uvec>(0,inputs_.n_cols-1,inputs_.n_cols);
        for (unsigned int i = 0; i < sample_size-1; ++i)
          std::swap (indices[i], indices[i + (rnd () % (inputs_.n_cols - i))]);
        for (unsigned int i = 1; i < sample_size; ++i)
          shuffled_indices_[i] = indices[i-1];
      }
      std::copy (shuffled_indices_.begin (), shuffled_indices_.begin () + sample_size, sample.begin ());
    }
    inline uint64_t
    rnd ()
    {
      return uint64_t(rnd_dist_(rnd_gen_));
    }
    virtual void getSamples(uint64_t iterations,arma::uvec&samples)
    {
        if (inputs_.n_cols < getSampleSize ())
        {
          std::cerr<<"Can not select"<<getSampleSize()<<"unique points out of"<<inputs_.n_cols;
          // one of these will make it stop :)
          samples.reset();
          iterations = INT_MAX - 1;
          return;
        }
        std::vector<int> samplesvec;
        samplesvec.resize(getSampleSize());
//        std::cerr<<"Get a second point which is different than the first"<<std::endl;
        for (unsigned int iter = 0; iter < max_sample_checks_; ++iter)
        {
//          std::cerr<<"Choose the random indices"<<std::endl;
          if (samples_radius_ < std::numeric_limits<double>::epsilon ())
          {
//              std::cerr<<"drawIndexSample"<<std::endl;
              drawIndexSample(samplesvec);
          }
          else
          {
//              std::cerr<<"drawIndexSampleRadius"<<std::endl;
              drawIndexSampleRadius(samplesvec);
          }

//          std::cerr<<"If it's a good sample, stop here"<<std::endl;
          if (isSampleGood (samplesvec))
          {
//            std::clog<<"Selected "<<samplesvec.size()<<" samples"<<std::endl;
            samples = arma::conv_to<arma::uvec>::from(samplesvec);
            return;
          }
        }
        std::cerr<<"Could not select"<<getSampleSize()<<"sample points in "<<max_sample_checks_<<" iterations"<<std::endl;
        samplesvec.clear ();
    }
    virtual bool computeModel(arma::uvec&inliers,arma::vec& coeff)=0;
    virtual void optimizeModel(const arma::uvec&inliers,
                               const arma::vec& coeff,
                               arma::vec& optimized_coeff)=0;
    virtual uint64_t countWithinDistance(arma::vec& coeff,double threshold)=0;
    virtual void selectWithinDistance(arma::vec& coeff,double threshold,arma::uvec& inliers)=0;
    virtual Model type()=0;
protected:
    virtual bool isSampleGood (const std::vector<int> &samples)const = 0;
    virtual inline bool isModelValid (const arma::vec& coeff) = 0;
    virtual int  getSampleSize() const=0;
protected:
    arma::fmat inputs_;
    arma::uvec shuffled_indices_;
    std::random_device rnd_d_;
    std::mt19937 rnd_gen_;
    std::uniform_int_distribution<> rnd_dist_;
    float samples_radius_;
    SearchPtr samples_radius_search_;
    static const unsigned int max_sample_checks_ = 1000;
};
}
#endif
