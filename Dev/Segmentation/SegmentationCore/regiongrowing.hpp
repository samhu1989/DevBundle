#include "regiongrowing.h"
#include "common.h"
#include <queue>
#include <list>
#include <cmath>
#include <time.h>
namespace Segmentation {
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template <typename M>
RegionGrowing<M>::RegionGrowing () :
  min_pts_per_cluster_ (1),
  max_pts_per_cluster_ (std::numeric_limits<int>::max ()),
  smooth_mode_flag_ (true),
  curvature_flag_ (true),
  residual_flag_ (false),
  theta_threshold_ (30.0f / 180.0f * static_cast<float> (M_PI)),
  residual_threshold_ (0.05f),
  curvature_threshold_ (0.05f),
  neighbour_number_ (30),
  neighbour_radius_ (0.0),
  search_ (),
  normals_(NULL),
  curvatures_(),
  input_(NULL),
  point_neighbours_ (0),
  point_labels_ (0),
  normal_flag_ (true),
  num_pts_in_segment_ (0),
  clusters_ (0),
  number_of_segments_ (0)
{
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template <typename M>
RegionGrowing<M>::~RegionGrowing ()
{
  if (search_ != 0)
    search_.reset ();
  if (normals_ != NULL)
    normals_=NULL;

  point_neighbours_.clear ();
  point_labels_.clear ();
  num_pts_in_segment_.clear ();
  clusters_.clear ();
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template <typename M> int
RegionGrowing<M>::getMinClusterSize ()
{
  return (min_pts_per_cluster_);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template <typename M> void
RegionGrowing<M>::setMinClusterSize (int min_cluster_size)
{
  min_pts_per_cluster_ = min_cluster_size;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template <typename M> int
RegionGrowing<M>::getMaxClusterSize ()
{
  return (max_pts_per_cluster_);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template <typename M> void
RegionGrowing<M>::setMaxClusterSize (int max_cluster_size)
{
  max_pts_per_cluster_ = max_cluster_size;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template <typename M> bool
RegionGrowing<M>::getSmoothModeFlag () const
{
  return (smooth_mode_flag_);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template <typename M> void
RegionGrowing<M>::setSmoothModeFlag (bool value)
{
  smooth_mode_flag_ = value;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template <typename M> bool
RegionGrowing<M>::getCurvatureTestFlag () const
{
  return (curvature_flag_);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template <typename M> void
RegionGrowing<M>::setCurvatureTestFlag (bool value)
{
  curvature_flag_ = value;

  if (curvature_flag_ == false && residual_flag_ == false)
    residual_flag_ = true;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template <typename M> bool
RegionGrowing<M>::getResidualTestFlag () const
{
  return (residual_flag_);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template <typename M> void
RegionGrowing<M>::setResidualTestFlag (bool value)
{
  residual_flag_ = value;

  if (curvature_flag_ == false && residual_flag_ == false)
    curvature_flag_ = true;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template <typename M> float
RegionGrowing<M>::getSmoothnessThreshold () const
{
  return (theta_threshold_);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template <typename M> void
RegionGrowing<M>::setSmoothnessThreshold (float theta)
{
  theta_threshold_ = theta;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template <typename M> float
RegionGrowing<M>::getResidualThreshold () const
{
  return (residual_threshold_);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template <typename M> void
RegionGrowing<M>::setResidualThreshold (float residual)
{
  residual_threshold_ = residual;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template <typename M> float
RegionGrowing<M>::getCurvatureThreshold () const
{
  return (curvature_threshold_);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template <typename M> void
RegionGrowing<M>::setCurvatureThreshold (float curvature)
{
  curvature_threshold_ = curvature;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template <typename M> unsigned int
RegionGrowing<M>::getNumberOfNeighbours () const
{
  return (neighbour_number_);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template <typename M> void
RegionGrowing<M>::setNumberOfNeighbours (unsigned int neighbour_number)
{
  neighbour_number_ = neighbour_number;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template <typename M> float
RegionGrowing<M>::getRadiusOfNeighbours () const
{
  return (neighbour_radius_);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template <typename M> void
RegionGrowing<M>::setRadiusOfNeighbours (float neighbour_radius)
{
  neighbour_radius_ = neighbour_radius;
  //when using neighbour_radius remove the limit of of neighbour_number
//  neighbour_number_ = std::numeric_limits<unsigned int>::max();
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template <typename M> typename RegionGrowing<M>::KdTreePtr
RegionGrowing<M>::getSearchMethod () const
{
  return (search_);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template <typename M> void
RegionGrowing<M>::setSearchMethod (const KdTreePtr& tree)
{
  if (search_ != 0)
    search_.reset ();

  search_ = tree;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template <typename M> typename RegionGrowing<M>::NormalPtr
RegionGrowing<M>::getInputNormals () const
{
  return (normals_);
}

template <typename M> typename RegionGrowing<M>::CurvaturePtr&
RegionGrowing<M>::getCurvatures()
{
    return curvatures_;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template <typename M> void
RegionGrowing<M>::setInputNormals (const NormalPtr& norm)
{
  if (normals_ != NULL)
    normals_ = NULL;

  normals_ = norm;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template <typename M> void
RegionGrowing<M>::extract (std::vector<arma::uvec>& clusters)
{
  clusters_.clear ();
  clusters.clear ();
  point_neighbours_.clear ();
  point_labels_.clear ();
  num_pts_in_segment_.clear ();
  number_of_segments_ = 0;

  bool segmentation_is_possible = initCompute ();
  if ( !segmentation_is_possible )
  {
    deinitCompute ();
    return;
  }

  segmentation_is_possible = prepareForSegmentation ();
  if ( !segmentation_is_possible )
  {
    deinitCompute ();
    return;
  }

  findPointNeighbours ();
  applySmoothRegionGrowingAlgorithm ();
  assembleRegions ();

  clusters.resize (clusters_.size ());
  std::vector<arma::uvec>::iterator cluster_iter_input = clusters.begin ();
  for (std::vector<arma::uvec>::const_iterator cluster_iter = clusters_.begin (); cluster_iter != clusters_.end (); cluster_iter++)
  {
    if ((static_cast<int> (cluster_iter->size ()) >= min_pts_per_cluster_) &&
        (static_cast<int> (cluster_iter->size ()) <= max_pts_per_cluster_))
    {
      *cluster_iter_input = *cluster_iter;
      cluster_iter_input++;
    }
  }

  clusters_ = std::vector<arma::uvec> (clusters.begin (), cluster_iter_input);
  clusters.resize(clusters_.size());

  deinitCompute ();
}
///////////////////////////////////////////////////////////////////////////////////
template <typename M>
void RegionGrowing<M>::extract (arma::uvec&labels)
{
    std::vector<arma::uvec> clusters;
    clusters_.clear ();
    clusters.clear ();
    point_neighbours_.clear ();
    point_labels_.clear ();
    num_pts_in_segment_.clear ();
    number_of_segments_ = 0;

    bool segmentation_is_possible = initCompute ();
    if ( !segmentation_is_possible )
    {
      deinitCompute ();
      return;
    }

    segmentation_is_possible = prepareForSegmentation ();
    if ( !segmentation_is_possible )
    {
      deinitCompute ();
      return;
    }

    findPointNeighbours ();
    applySmoothRegionGrowingAlgorithm ();
    assembleRegions ();

    clusters.resize (clusters_.size ());
    std::vector<arma::uvec>::iterator cluster_iter_input = clusters.begin ();
    for (std::vector<arma::uvec>::const_iterator cluster_iter = clusters_.begin (); cluster_iter != clusters_.end (); cluster_iter++)
    {
      if ((static_cast<int> (cluster_iter->size ()) >= min_pts_per_cluster_) &&
          (static_cast<int> (cluster_iter->size ()) <= max_pts_per_cluster_))
      {
        *cluster_iter_input = *cluster_iter;
        cluster_iter_input++;
      }
    }

    clusters_ = std::vector<arma::uvec> (clusters.begin (), cluster_iter_input);
    clusters.resize(clusters_.size());

    arma::uword baselabel=0;
    if(indices_.size()!=labels.size())
    {
        arma::uvec tmplabel = labels;
        tmplabel.elem(indices_).fill(0);
        baselabel = arma::max(tmplabel);
    }

    int l;
    for(l=0;l<clusters.size();++l)
    {
        labels.elem(clusters[l]).fill(arma::uword(baselabel+l+1));
    }

    deinitCompute ();
}
///////////////////////////////////////////////////////////////////////////////////
template <typename M>
void RegionGrowing<M>::setInputMesh(MeshPtr input)
{
    input_ = static_cast<MeshPtr>(input);
    if(input_)
    {
        if(input_->has_vertex_normals())
        {
            normals_=(float*)input_->vertex_normals();
        }
        if(indices_.is_empty())search_cloud_.reset(new MeshKDTreeInterface<M>(*input_));
        else search_cloud_.reset(new MeshKDTreeInterface<M>(*input_,indices_));
    }else{
        std::logic_error("Input NULL");
    }
}

template <typename M>
void RegionGrowing<M>::setIndices(arma::uvec&indices)
{
    indices_ = indices;
    if(input_)search_cloud_.reset(new MeshKDTreeInterface<M>(*input_,indices_));
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template <typename M> bool
RegionGrowing<M>::prepareForSegmentation ()
{
  // if user forgot to pass point cloud or if it is empty
  if(!input_)return false;

  if ( input_->n_vertices() == 0 )
    return (false);

  if( indices_.size() == 0)
  {
      indices_ = arma::linspace<arma::uvec>(0,input_->n_vertices()-1,input_->n_vertices());
  }

  // if user forgot to pass normals or the sizes of point and normal cloud are different
  if ( !normals_ )
    return (false);

  // if residual test is on then we need to check if all needed parameters were correctly initialized
  if (residual_flag_)
  {
    if (residual_threshold_ <= 0.0f)
      return (false);
  }

  // if curvature test is on ...
   if (curvature_flag_||normal_flag_)
   {
     if(!curvatures_)return false;
   }

  // from here we check those parameters that are always valuable
  if (neighbour_number_ == 0)
    return (false);

  if( neighbour_number_ >= input_->n_vertices())
      return false;

  // if user didn't set search method
  if (!search_)
  {
    search_.reset (new KdTree(3,*search_cloud_,nanoflann::KDTreeSingleIndexAdaptorParams(3)));
    search_->buildIndex();
  }

  return (true);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template <typename M> void
RegionGrowing<M>::findPointNeighbours ()
{
  int point_number = static_cast<int> (indices_.size());
  point_neighbours_.resize (input_->n_vertices());
  float* pts_ptr = (float*)input_->points();
  std::vector<std::pair<arma::uword,float>> radiusResult;
  arma::uword* knnNeighbors = new arma::uword[neighbour_number_];
  float* knnDistance = new float[neighbour_number_];
  nanoflann::SearchParams param;
  for (int i = 0; i < point_number; i++)
  {
      int pi = indices_(i);
      if (!std::isfinite(pts_ptr[3*pi])||!std::isfinite(pts_ptr[3*pi+1])||!std::isfinite(pts_ptr[3*pi+2]))
          continue;
      if(neighbour_number_==0){
          search_->radiusSearch(&pts_ptr[3*pi],neighbour_radius_,radiusResult,param);
          std::vector<std::pair<arma::uword,float>>::iterator iter;
          point_neighbours_[pi].clear();
          for(iter=radiusResult.begin();iter!=radiusResult.end();++iter)
          {
              point_neighbours_[pi].push_back(indices_(iter->first));
          }
      }
      else {
          search_->knnSearch(&pts_ptr[3*pi], neighbour_number_,knnNeighbors,knnDistance);
          for(int k=0;k<neighbour_number_;++k)
          {
              if( std::sqrt(knnDistance[k]) <= neighbour_radius_ )point_neighbours_[pi].push_back(indices_(knnNeighbors[k]));
          }
      }
  }
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template <typename M> void
RegionGrowing<M>::applySmoothRegionGrowingAlgorithm ()
{
  int num_of_pts = static_cast<int> (indices_.size());
  point_labels_.resize (input_->n_vertices(), -1);

  std::vector< std::pair<float, int> > point_residual;
  std::pair<float, int> pair;
  point_residual.resize (num_of_pts, pair);

  float* c_ptr = curvatures_.get();

  if (normal_flag_ == true)
  {
    for (int i_point = 0; i_point < num_of_pts; i_point++)
    {
      int pi = indices_(i_point);
      point_residual[i_point].first = c_ptr[pi];
      point_residual[i_point].second = pi;
    }
    std::sort (point_residual.begin (), point_residual.end (), comparePair);
  }
  else
  {
    for (int i_point = 0; i_point < num_of_pts; i_point++)
    {
      int pi = indices_(i_point);
      point_residual[i_point].first = 0;
      point_residual[i_point].second = pi;
    }
  }
  int seed_counter = 0;
  int seed = point_residual[seed_counter].second;

  int segmented_pts_num = 0;
  int number_of_segments = 0;
  while (segmented_pts_num < num_of_pts)
  {
    int pts_in_segment;
    pts_in_segment = growRegion (seed, number_of_segments);
    segmented_pts_num += pts_in_segment;
    num_pts_in_segment_.push_back (pts_in_segment);
    number_of_segments++;

    //find next point that is not segmented yet
    for (int i_seed = seed_counter + 1; i_seed < num_of_pts; i_seed++)
    {
      int index = point_residual[i_seed].second;
      if (point_labels_[index] == -1)
      {
        seed = index;
        break;
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template <typename M> int
RegionGrowing<M>::growRegion (int initial_seed, int segment_number)
{
  std::queue<int> seeds;
  seeds.push (initial_seed);
  point_labels_[initial_seed] = segment_number;

  int num_pts_in_segment = 1;

  while (!seeds.empty ())
  {
    int curr_seed;
    curr_seed = seeds.front ();
    seeds.pop ();

    size_t i_nghbr = 0;
    while ( i_nghbr < neighbour_number_ && i_nghbr < point_neighbours_[curr_seed].size () )
    {
      int index = point_neighbours_[curr_seed][i_nghbr];
      if (point_labels_[index] != -1)
      {
        i_nghbr++;
        continue;
      }

      bool is_a_seed = false;
      bool belongs_to_segment = validatePoint (initial_seed, curr_seed, index, is_a_seed);

      if (belongs_to_segment == false)
      {
        i_nghbr++;
        continue;
      }

      point_labels_[index] = segment_number;
      num_pts_in_segment++;

      if (is_a_seed)
      {
        seeds.push (index);
      }

      i_nghbr++;
    }// next neighbour
  }// next seed

  return (num_pts_in_segment);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template <typename M> bool
RegionGrowing<M>::validatePoint (int initial_seed, int point, int nghbr, bool& is_a_seed) const
{
  is_a_seed = true;

  float cosine_threshold = cosf (theta_threshold_);
  float* pts_ptr = (float*)input_->points();
  float* n_ptr = normals_;

  arma::fvec initial_point (&(pts_ptr[3*point]),3,false,true);
  arma::fvec initial_normal (&(n_ptr[3*point]),3,false,true);

  //check the angle between normals
  if (smooth_mode_flag_ == true)
  {
    arma::fvec nghbr_normal (&(n_ptr[3*nghbr]),3,false,true);
    float dot_product = std::fabs(arma::dot(initial_normal,nghbr_normal));
    if (dot_product < cosine_threshold)
    {
      return (false);
    }
  }
  else
  {
    arma::fvec nghbr_normal (&(n_ptr[3*nghbr]),3,false,true);
    arma::fvec initial_seed_normal(&(n_ptr[initial_seed]),3,false,true);
    float dot_product = std::fabs(arma::dot(initial_seed_normal,nghbr_normal));
    if (dot_product < cosine_threshold)
      return (false);
  }

  // check the curvature if needed
  float* c_ptr = curvatures_.get();
  if (curvature_flag_ && c_ptr[nghbr] > curvature_threshold_)
  {
    is_a_seed = false;
  }

  // check the residual if needed
  arma::fvec nghbr_point(&pts_ptr[3*nghbr],3,false,true);
  float residual = std::fabs(arma::dot(initial_normal,initial_point - nghbr_point));
  if (residual_flag_ && residual > residual_threshold_)
    is_a_seed = false;

  return (true);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template <typename M> void
RegionGrowing<M>::assembleRegions ()
{
  int number_of_segments = static_cast<int> (num_pts_in_segment_.size ());
  int number_of_points = static_cast<int> (input_->n_vertices());

  arma::uvec segment;
  clusters_.resize (number_of_segments, segment);

  for (int i_seg = 0; i_seg < number_of_segments; i_seg++)
  {
    clusters_[i_seg] = arma::uvec(num_pts_in_segment_[i_seg],arma::fill::zeros);
  }

  std::vector<int> counter;
  counter.resize (number_of_segments, 0);

  for (int i_point = 0; i_point < number_of_points; i_point++)
  {
    int segment_index = point_labels_[i_point];
    if (segment_index != -1)
    {
      int point_index = counter[segment_index];
      clusters_[segment_index][point_index] = i_point;
      counter[segment_index] = point_index + 1;
    }
  }

  number_of_segments_ = number_of_segments;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template <typename M> void
RegionGrowing<M>::getSegmentFromPoint (int index, arma::uvec &cluster)
{
  cluster.reset();

  bool segmentation_is_possible = initCompute ();
  if ( !segmentation_is_possible )
  {
    deinitCompute ();
    return;
  }

  // first of all we need to find out if this point belongs to cloud
  bool point_was_found = false;
  int number_of_points = static_cast <int> (indices_.size());
  for (int point = 0; point < number_of_points; point++)
    if ( indices_(point) == index)
    {
      point_was_found = true;
      break;
    }

  if (point_was_found)
  {
    if (clusters_.empty ())
    {
      point_neighbours_.clear ();
      point_labels_.clear ();
      num_pts_in_segment_.clear ();
      number_of_segments_ = 0;

      segmentation_is_possible = prepareForSegmentation ();
      if ( !segmentation_is_possible )
      {
        deinitCompute ();
        return;
      }

      findPointNeighbours ();
      applySmoothRegionGrowingAlgorithm ();
      assembleRegions ();
    }
    // if we have already made the segmentation, then find the segment
    // to which this point belongs
    std::vector <arma::uvec>::iterator i_segment;
    for (i_segment = clusters_.begin (); i_segment != clusters_.end (); i_segment++)
    {
      bool segment_was_found = false;
      for (size_t i_point = 0; i_point < i_segment->size (); i_point++)
      {
        if ((*i_segment)[i_point] == index)
        {
          segment_was_found = true;
          cluster = *i_segment;
          break;
        }
      }
      if (segment_was_found)
      {
        break;
      }
    }// next segment
  }// end if point was found

  deinitCompute ();
}

template <typename M>
bool RegionGrowing<M>::initCompute()
{
    return true;
}

template <typename M>
void RegionGrowing<M>::deinitCompute()
{
    curvatures_.reset();
    normals_ = NULL;
    search_.reset();
    search_cloud_.reset();
    point_neighbours_.clear ();
    point_labels_.clear ();
    num_pts_in_segment_.clear ();
    clusters_.clear ();
    input_ = NULL;
    indices_.reset();
}

}
