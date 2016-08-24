#ifndef REGIONGROWINGRGB_HPP
#define REGIONGROWINGRGB_HPP
#include "regiongrowingrgb.h"
#include <queue>
#include <cmath>

namespace Segmentation{
template<typename M>
RegionGrowingRGB<M>::RegionGrowingRGB():
    RegionGrowing<M>(),
    color_p2p_threshold_ (1225.0f),
    color_r2r_threshold_ (10.0f),
    distance_threshold_ (0.05f),
    region_neighbour_number_ (100),
    point_distances_ (0),
    segment_neighbours_ (0),
    segment_distances_ (0),
    segment_labels_ (0)
{
    normal_flag_ = false;
    curvature_flag_ = false;
    residual_flag_ = false;
    min_pts_per_cluster_ = 10;
}

template<typename M>
RegionGrowingRGB<M>::~RegionGrowingRGB()
{
    point_distances_.clear ();
    segment_neighbours_.clear ();
    segment_distances_.clear ();
    segment_labels_.clear ();
}

template <typename M> float
RegionGrowingRGB<M>::getPointColorThreshold () const
{
    return (std::pow(color_p2p_threshold_, 0.5f));
}

template <typename M> void
RegionGrowingRGB<M>::setPointColorThreshold(float thresh)
{
    color_p2p_threshold_ = thresh * thresh;
}

template <typename M> void
RegionGrowingRGB<M>::setRegionColorThreshold (float thresh)
{
    color_r2r_threshold_ = thresh * thresh;
}

template <typename M> float
RegionGrowingRGB<M>::getRegionColorThreshold () const
{
    return (std::pow(color_r2r_threshold_, 0.5f));
}

template <typename M> float
RegionGrowingRGB<M>::getDistanceThreshold () const
{
    return (std::pow(distance_threshold_, 0.5f));
}

template <typename M> void
RegionGrowingRGB<M>::setDistanceThreshold (float thresh)
{
    distance_threshold_ = thresh * thresh;
}

template <typename M> unsigned int
RegionGrowingRGB<M>::getNumberOfRegionNeighbours () const
{
    return (region_neighbour_number_);
}

template <typename M> void
RegionGrowingRGB<M>::setNumberOfRegionNeighbours(unsigned int nghbr_number)
{
    region_neighbour_number_ = nghbr_number;
}

template <typename M> bool
RegionGrowingRGB<M>::getNormalTestFlag () const
{
    return (normal_flag_);
}

template <typename M> void
RegionGrowingRGB<M>::setNormalTestFlag (bool value)
{
    normal_flag_ = value;
}

template <typename M> void
RegionGrowingRGB<M>::setCurvatureTestFlag (bool value)
{
    curvature_flag_ = value;
}

template <typename M> void
RegionGrowingRGB<M>::setResidualTestFlag (bool value)
{
    residual_flag_ = value;
}

template <typename M> void
RegionGrowingRGB<M>::extract (std::vector<arma::uvec>& clusters)
{
  clusters_.clear ();
  clusters.clear ();
  point_neighbours_.clear ();
  point_labels_.clear ();
  num_pts_in_segment_.clear ();
  point_distances_.clear ();
  segment_neighbours_.clear ();
  segment_distances_.clear ();
  segment_labels_.clear ();
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

  findPointNeighbours();
  applySmoothRegionGrowingAlgorithm();
  RegionGrowing<M>::assembleRegions();

  findSegmentNeighbours();
  applyRegionMergingAlgorithm();

  std::vector<arma::uvec>::iterator cluster_iter = clusters_.begin ();
  while (cluster_iter != clusters_.end ())
  {
       if (static_cast<int> (cluster_iter->size()) < min_pts_per_cluster_ ||
           static_cast<int> (cluster_iter->size()) > max_pts_per_cluster_)
       {
            cluster_iter = clusters_.erase (cluster_iter);
       }
       else
            cluster_iter++;
  }

  clusters.reserve (clusters_.size ());
  std::copy (clusters_.begin (), clusters_.end (), std::back_inserter (clusters));
  deinitCompute ();
}

template<typename M> void
RegionGrowingRGB<M>::extract(arma::uvec& labels)
{
  std::vector<arma::uvec> clusters;
  clusters_.clear ();
  clusters.clear();
  point_neighbours_.clear ();
  point_labels_.clear ();
  num_pts_in_segment_.clear ();
  point_distances_.clear ();
  segment_neighbours_.clear ();
  segment_distances_.clear ();
  segment_labels_.clear ();
  number_of_segments_ = 0;

//  std::cerr<<"0"<<std::endl;

  bool segmentation_is_possible = initCompute ();
  if ( !segmentation_is_possible )
  {
    deinitCompute ();
    return;
  }

//  std::cerr<<"1"<<std::endl;

  segmentation_is_possible = prepareForSegmentation ();
  if ( !segmentation_is_possible )
  {
    deinitCompute ();
    return;
  }

//  std::cerr<<"2"<<std::endl;

  findPointNeighbours();
  applySmoothRegionGrowingAlgorithm();
  RegionGrowing<M>::assembleRegions();


//  std::cerr<<"clusters_.size():"<<clusters_.size()<<std::endl;
  findSegmentNeighbours();
  applyRegionMergingAlgorithm();

//  std::cerr<<"4"<<std::endl;

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
//  std::cerr<<"6"<<std::endl;
  deinitCompute ();
}

template <typename M> bool
RegionGrowingRGB<M>::prepareForSegmentation ()
{
  // if user forgot to pass point cloud or if it is empty
  if ( input_->n_vertices() == 0 )
  return (false);

  // if normal/smoothness test is on then we need to check if all needed variables and parameters
  // were correctly initialized
  if (normal_flag_)
  {
    // if user forgot to pass normals or the sizes of point and normal cloud are different
    if ( normals_ == NULL )
    return (false);
  }

  // if residual test is on then we need to check if all needed parameters were correctly initialized
  if (residual_flag_)
  {
    if (residual_threshold_ <= 0.0f)
    return (false);
  }

  // if curvature test is on ...
//  if (curvature_flag_||normal_flag_)
//  {
//    if(!curvatures_)return false;
//  }

  // here we check the parameters related to color-based segmentation
  if ( region_neighbour_number_ == 0 || color_p2p_threshold_ < 0.0f || color_r2r_threshold_ < 0.0f || distance_threshold_ < 0.0f )
  return (false);

  // from here we check those parameters that are always valuable
  if (neighbour_number_ == 0)
  return (false);

  // if user didn't set search method
  if (!search_)
  {
    search_.reset (new KdTree(3,*search_cloud_,nanoflann::KDTreeSingleIndexAdaptorParams(3)));
    search_->buildIndex();
  }

  if( indices_.size() == 0)
  {
      indices_ = arma::linspace<arma::uvec>(0,input_->n_vertices()-1,input_->n_vertices());
  }
  return (true);
}

template <typename M> void
RegionGrowingRGB<M>::findPointNeighbours ()
{
    int point_number = static_cast<int> (indices_.size());
    point_neighbours_.resize (input_->n_vertices());
    point_distances_.resize (input_->n_vertices());
    float* pts_ptr = (float*)input_->points();
//    std::vector<std::pair<arma::uword,float>> radiusResult;
    arma::uword* knnNeighbors = new arma::uword[neighbour_number_];
    float* knnDistance = new float[neighbour_number_];
    nanoflann::SearchParams param;
    for (int i = 0; i < point_number; i++)
    {
        int pi = indices_(i);
        if (!std::isfinite(pts_ptr[3*pi])||!std::isfinite(pts_ptr[3*pi+1])||!std::isfinite(pts_ptr[3*pi+2]))
            continue;
        assert(neighbour_number_>0);
        search_->knnSearch(&pts_ptr[3*pi], neighbour_number_,knnNeighbors,knnDistance);
        for(int k=0;k<neighbour_number_;++k)
        {
            point_neighbours_[pi].push_back(indices_(knnNeighbors[k]));
            point_distances_[pi].push_back(knnDistance[k]);
        }
    }
}

template <typename M> void
RegionGrowingRGB<M>::findSegmentNeighbours()
{
//  std::cerr<<"RegionGrowingRGB<M>::findSegmentNeighbours(0)"<<std::endl;
//  std::cerr<<"number of segments:"<<number_of_segments_<<std::endl;
  std::vector<int> neighbours;
  std::vector<float> distances;
  segment_neighbours_.resize (number_of_segments_, neighbours);
  segment_distances_.resize (number_of_segments_, distances);
//  std::cerr<<"RegionGrowingRGB<M>::findSegmentNeighbours(1)"<<std::endl;
  for (int i_seg = 0; i_seg < number_of_segments_; i_seg++)
  {
       std::vector<int> nghbrs;
       std::vector<float> dist;
       findRegionsKNN (i_seg, region_neighbour_number_, nghbrs, dist);
       segment_neighbours_[i_seg].swap (nghbrs);
       segment_distances_[i_seg].swap (dist);
  }
//  std::cerr<<"RegionGrowingRGB<M>::findSegmentNeighbours(2)"<<std::endl;
}

template<typename M> void
RegionGrowingRGB<M>::findRegionsKNN(
        int index,
        int nghbr_number,
        std::vector<int>& nghbrs,
        std::vector<float>& dist
        )
{
//  std::cerr<<"RegionGrowingRGB<M>::findRegionsKNN(0)"<<std::endl;
  std::vector<float> distances;
  float max_dist = std::numeric_limits<float>::max ();
  distances.resize (clusters_.size (), max_dist);

  int number_of_points = num_pts_in_segment_[index];
  //loop throug every point in this segment and check neighbours
//  std::cerr<<"num_pts_in_segment_.size():"<<num_pts_in_segment_.size()<<std::endl;
//  std::cerr<<"clusters_.size():"<<clusters_.size()<<std::endl;
//  std::cerr<<"index:"<<index<<std::endl;
//  std::cerr<<"point_labels_:"<<point_labels_.size()<<std::endl;
  for (int i_point = 0; i_point < number_of_points; i_point++)
  {
    assert( i_point < clusters_[index].size());
    int point_index = static_cast<int>(clusters_[index](i_point));
    assert( point_index < point_neighbours_.size() );
    int number_of_neighbours = static_cast<int>(point_neighbours_[point_index].size ());
    //loop throug every neighbour of the current point, find out to which segment it belongs
    //and if it belongs to neighbouring segment and is close enough then remember segment and its distance
//    std::cerr<<"RegionGrowingRGB<M>::findRegionsKNN(0.1)"<<std::endl;
    for (int i_nghbr = 0; i_nghbr < number_of_neighbours; i_nghbr++)
    {
        // find segment
        int segment_index = -1;
        segment_index = point_labels_[point_neighbours_[point_index][i_nghbr]];
        assert( (segment_index < distances.size()) );
        if ( (segment_index != index) && ( segment_index > 0 ) )
        {
           // try to push it to the queue
           assert( point_index < point_distances_.size() );
           assert( i_nghbr < point_distances_[point_index].size() );
           if ( distances[segment_index] > point_distances_[point_index][i_nghbr] )
             distances[segment_index] = point_distances_[point_index][i_nghbr];
        }
     }
//    std::cerr<<"RegionGrowingRGB<M>::findRegionsKNN(0.2)"<<std::endl;
  }// next point
//  std::cerr<<"RegionGrowingRGB<M>::findRegionsKNN(1)"<<std::endl;

  std::priority_queue<std::pair<float, int> > segment_neighbours;
  for (int i_seg = 0; i_seg < number_of_segments_; i_seg++)
  {
    if (distances[i_seg] < max_dist)
    {
        segment_neighbours.push (std::make_pair (distances[i_seg], i_seg) );
        if (int (segment_neighbours.size ()) > nghbr_number)
            segment_neighbours.pop ();
    }
  }
//  std::cerr<<"RegionGrowingRGB<M>::findRegionsKNN(2)"<<std::endl;
  int size = std::min<int> (static_cast<int> (segment_neighbours.size ()), nghbr_number);
  nghbrs.resize (size, 0);
  dist.resize (size, 0);
  int counter = 0;
  while ( !segment_neighbours.empty () && counter < nghbr_number )
  {
    dist[counter] = segment_neighbours.top().first;
    nghbrs[counter] = segment_neighbours.top().second;
    segment_neighbours.pop();
    counter++;
  }
//  std::cerr<<"RegionGrowingRGB<M>::findRegionsKNN(3)"<<std::endl;
}

template <typename M> void
RegionGrowingRGB<M>::applyRegionMergingAlgorithm ()
{
  int number_of_points = static_cast<int> (indices_.size ());

  // calculate color of each segment
  std::vector< std::vector<unsigned int> > segment_color;
  std::vector<unsigned int> color;
  color.resize (3, 0);
  segment_color.resize (number_of_segments_, color);
  uint8_t* colors = (uint8_t*)input_->vertex_colors();
  for (int i_point = 0; i_point < number_of_points; i_point++)
  {
       int point_index = indices_(i_point);
       int segment_index = point_labels_[point_index];
       segment_color[segment_index][0] += (unsigned int)(colors[3*point_index+0]);
       segment_color[segment_index][1] += (unsigned int)(colors[3*point_index+1]);
       segment_color[segment_index][2] += (unsigned int)(colors[3*point_index+2]);
  }
  for (int i_seg = 0; i_seg < number_of_segments_; i_seg++)
  {
    segment_color[i_seg][0] = static_cast<unsigned int> (static_cast<float> (segment_color[i_seg][0]) / static_cast<float> (num_pts_in_segment_[i_seg]));
    segment_color[i_seg][1] = static_cast<unsigned int> (static_cast<float> (segment_color[i_seg][1]) / static_cast<float> (num_pts_in_segment_[i_seg]));
    segment_color[i_seg][2] = static_cast<unsigned int> (static_cast<float> (segment_color[i_seg][2]) / static_cast<float> (num_pts_in_segment_[i_seg]));
  }

  // now it is time to find out if there are segments with a similar color
  // and merge them together
  std::vector<unsigned int> num_pts_in_homogeneous_region;
  std::vector<int> num_seg_in_homogeneous_region;

  segment_labels_.resize (number_of_segments_, -1);

  float dist_thresh = distance_threshold_;
  int homogeneous_region_number = 0;
  int curr_homogeneous_region = 0;
  for (int i_seg = 0; i_seg < number_of_segments_; i_seg++)
  {
      curr_homogeneous_region = 0;
      if (segment_labels_[i_seg] == -1)
      {
         segment_labels_[i_seg] = homogeneous_region_number;
         curr_homogeneous_region = homogeneous_region_number;
         num_pts_in_homogeneous_region.push_back (num_pts_in_segment_[i_seg]);
         num_seg_in_homogeneous_region.push_back (1);
         homogeneous_region_number++;
       }
       else
         curr_homogeneous_region = segment_labels_[i_seg];

       unsigned int i_nghbr = 0;
       while ( i_nghbr < region_neighbour_number_ && i_nghbr < segment_neighbours_[i_seg].size () )
       {
         int index = segment_neighbours_[i_seg][i_nghbr];
         if (segment_distances_[i_seg][i_nghbr] > dist_thresh)
         {
           i_nghbr++;
           continue;
         }
         if ( segment_labels_[index] == -1 )
         {
           float difference = calculateColorimetricalDifference (segment_color[i_seg], segment_color[index]);
           if (difference < color_r2r_threshold_)
           {
             segment_labels_[index] = curr_homogeneous_region;
             num_pts_in_homogeneous_region[curr_homogeneous_region] += num_pts_in_segment_[index];
             num_seg_in_homogeneous_region[curr_homogeneous_region] += 1;
           }
         }
         i_nghbr++;
       }// next neighbour
     }// next segment

     segment_color.clear ();
     color.clear ();

     std::vector< std::vector<int> > final_segments;
     std::vector<int> region;
     final_segments.resize (homogeneous_region_number, region);
     for (int i_reg = 0; i_reg < homogeneous_region_number; i_reg++)
     {
       final_segments[i_reg].resize (num_seg_in_homogeneous_region[i_reg], 0);
     }

  std::vector<int> counter;
  counter.resize (homogeneous_region_number, 0);
  for (int i_seg = 0; i_seg < number_of_segments_; i_seg++)
  {
    int index = segment_labels_[i_seg];
    final_segments[ index ][ counter[index] ] = i_seg;
    counter[index] += 1;
  }

  std::vector< std::vector< std::pair<float, int> > > region_neighbours;
  findRegionNeighbours (region_neighbours, final_segments);

        int final_segment_number = homogeneous_region_number;
        for (int i_reg = 0; i_reg < homogeneous_region_number; i_reg++)
        {
            if (static_cast<int> (num_pts_in_homogeneous_region[i_reg]) < min_pts_per_cluster_)
            {
                if ( region_neighbours[i_reg].empty () )
                continue;
         int nearest_neighbour = region_neighbours[i_reg][0].second;
         if ( region_neighbours[i_reg][0].first == std::numeric_limits<float>::max () )
           continue;
         int reg_index = segment_labels_[nearest_neighbour];
         int num_seg_in_reg = num_seg_in_homogeneous_region[i_reg];
         for (int i_seg = 0; i_seg < num_seg_in_reg; i_seg++)
         {
           int segment_index = final_segments[i_reg][i_seg];
           final_segments[reg_index].push_back (segment_index);
           segment_labels_[segment_index] = reg_index;
         }
         final_segments[i_reg].clear ();
         num_pts_in_homogeneous_region[reg_index] += num_pts_in_homogeneous_region[i_reg];
         num_pts_in_homogeneous_region[i_reg] = 0;
         num_seg_in_homogeneous_region[reg_index] += num_seg_in_homogeneous_region[i_reg];
         num_seg_in_homogeneous_region[i_reg] = 0;
         final_segment_number -= 1;

         int nghbr_number = static_cast<int> (region_neighbours[reg_index].size ());
         for (int i_nghbr = 0; i_nghbr < nghbr_number; i_nghbr++)
         {
           if ( segment_labels_[ region_neighbours[reg_index][i_nghbr].second ] == reg_index )
           {
             region_neighbours[reg_index][i_nghbr].first = std::numeric_limits<float>::max ();
             region_neighbours[reg_index][i_nghbr].second = 0;
           }
         }
         nghbr_number = static_cast<int> (region_neighbours[i_reg].size ());
         for (int i_nghbr = 0; i_nghbr < nghbr_number; i_nghbr++)
         {
           if ( segment_labels_[ region_neighbours[i_reg][i_nghbr].second ] != reg_index )
           {
             std::pair<float, int> pair;
             pair.first = region_neighbours[i_reg][i_nghbr].first;
             pair.second = region_neighbours[i_reg][i_nghbr].second;
             region_neighbours[reg_index].push_back (pair);
           }
         }
            region_neighbours[i_reg].clear ();
            std::sort (region_neighbours[reg_index].begin (), region_neighbours[reg_index].end (), comparePair);
        }
    }
    assembleRegions (num_pts_in_homogeneous_region, static_cast<int> (num_pts_in_homogeneous_region.size ()));
    number_of_segments_ = final_segment_number;
}

template <typename M> float
RegionGrowingRGB<M>::calculateColorimetricalDifference (std::vector<unsigned int>& first_color, std::vector<unsigned int>& second_color) const
{
  float difference = 0.0f;
  arma::Col<uint8_t> first_rgb = arma::conv_to<arma::Col<uint8_t>>::from(first_color);
  arma::Col<uint8_t> second_rgb = arma::conv_to<arma::Col<uint8_t>>::from(second_color);
  arma::fvec::fixed<3> first_Lab,second_Lab;
  ColorArray::RGB2Lab(first_rgb,first_Lab);
  ColorArray::RGB2Lab(second_rgb,second_Lab);
  difference += float ( ( first_Lab(1) - second_Lab(1) ) * (first_Lab(1) - second_Lab(1)));
  difference += float ((first_Lab(2) - second_Lab(2)) * (first_Lab(2) - second_Lab(2)));
  return (difference);
}

template<typename M> void
RegionGrowingRGB<M>::findRegionNeighbours(std::vector<std::vector<std::pair<float, int>>>& neighbours_out, std::vector<std::vector<int>>& regions_in)
{
    int region_number = static_cast<int> (regions_in.size ());
    neighbours_out.clear ();
    neighbours_out.resize (region_number);

    for (int i_reg = 0; i_reg < region_number; i_reg++)
    {
     int segment_num = static_cast<int> (regions_in[i_reg].size ());
     neighbours_out[i_reg].reserve (segment_num * region_neighbour_number_);
   for (int i_seg = 0; i_seg < segment_num; i_seg++)
     {
       int curr_segment = regions_in[i_reg][i_seg];
       int nghbr_number = static_cast<int> (segment_neighbours_[curr_segment].size ());
       std::pair<float, int> pair;
       for (int i_nghbr = 0; i_nghbr < nghbr_number; i_nghbr++)
       {
         int segment_index = segment_neighbours_[curr_segment][i_nghbr];
         if ( segment_distances_[curr_segment][i_nghbr] == std::numeric_limits<float>::max () )
           continue;
         if (segment_labels_[segment_index] != i_reg)
         {
           pair.first = segment_distances_[curr_segment][i_nghbr];
           pair.second = segment_index;
           neighbours_out[i_reg].push_back(pair);
         }
       }// next neighbour
     }// next segment
     std::sort (neighbours_out[i_reg].begin(),neighbours_out[i_reg].end(),comparePair);
   }// next homogeneous region
 }

 //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template <typename M> void
RegionGrowingRGB<M>::assembleRegions(std::vector<unsigned int>& num_pts_in_region, int num_regions)
{
    clusters_.clear ();
    arma::uvec segment;
    clusters_.resize (num_regions, segment);
    for (int i_seg = 0; i_seg < num_regions; i_seg++)
    {
        clusters_[i_seg].resize (num_pts_in_region[i_seg]);
    }

    std::vector<int> counter;
    counter.resize (num_regions, 0);
    int point_number = static_cast<int> (indices_.size ());
    for (int i_point = 0; i_point < point_number; i_point++)
    {
         int point_index = indices_(i_point);
         int index = point_labels_[point_index];
         index = segment_labels_[index];
         clusters_[index][ counter[index]] = point_index;
         counter[index] += 1;
    }

   // now we need to erase empty regions
   if (clusters_.empty ())
     return;

   std::vector<arma::uvec>::iterator itr1, itr2;
   itr1 = clusters_.begin ();
   itr2 = clusters_.end () - 1;

   while (itr1 < itr2)
   {
     while (!(itr1->empty ()) && itr1 < itr2)
       itr1++;
     while (  itr2->empty ()  && itr1 < itr2)
       itr2--;

     if (itr1 != itr2)
       itr1->swap(*itr2);
    }

     if (itr2->empty ())
     clusters_.erase (itr2, clusters_.end ());
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template <typename M> bool
RegionGrowingRGB<M>::validatePoint (int initial_seed, int point, int nghbr, bool& is_a_seed) const
{
    is_a_seed = true;
    // check the color difference
    std::vector<unsigned int> point_color;
    point_color.resize (3, 0);
    std::vector<unsigned int> nghbr_color;
    nghbr_color.resize (3, 0);
   uint8_t* color = (uint8_t*)input_->vertex_colors();
   point_color[0] = (unsigned int)(color[3*point+0]);
   point_color[1] = (unsigned int)(color[3*point+1]);
   point_color[2] = (unsigned int)(color[3*point+2]);
   nghbr_color[0] = (unsigned int)(color[3*nghbr+0]);
   nghbr_color[1] = (unsigned int)(color[3*nghbr+1]);
   nghbr_color[2] = (unsigned int)(color[3*nghbr+2]);
   float difference = calculateColorimetricalDifference (point_color, nghbr_color);
   if (difference > color_p2p_threshold_)
     return (false);

   float cosine_threshold = cosf (theta_threshold_);

   // check the angle between normals if needed
   if (normal_flag_)
   {
     float data[3];
     float *p = (float*)input_->points();
     data[0] = p[3*point+0];
     data[1] = p[3*point+1];
     data[2] = p[3*point+2];

     arma::fvec::fixed<3> initial_point (static_cast<float*> (data));
     arma::fvec::fixed<3> initial_normal (static_cast<float*>(normals_+3*point));
     if (smooth_mode_flag_ == true)
    {
       arma::fvec::fixed<3> nghbr_normal (static_cast<float*>(normals_+3*nghbr));
       float dot_product = fabs(arma::dot(nghbr_normal,initial_normal));
       if (dot_product < cosine_threshold)
         return (false);
     }
     else
     {
       arma::fvec::fixed<3> nghbr_normal (static_cast<float*> (normals_+3*nghbr));
       arma::fvec::fixed<3> initial_seed_normal (static_cast<float*> (normals_+3*initial_seed));
       float dot_product = fabs(arma::dot(nghbr_normal,initial_seed_normal));
       if (dot_product < cosine_threshold)
         return (false);
     }
   }

   // check the curvature if needed
   if (curvature_flag_ && curvatures_.get()[3*nghbr] > curvature_threshold_)
     is_a_seed = false;

   // check the residual if needed
   if (residual_flag_)
   {
     float *p = (float*)input_->points();
     float data_p[3];
     data_p[0] = p[3*point+0];
     data_p[1] = p[3*point+1];
     data_p[2] = p[3*point+2];
     float data_n[3];
     data_n[0] = p[3*nghbr+0];
     data_n[1] = p[3*nghbr+1];
     data_n[2] = p[3*nghbr+2];
     arma::fvec::fixed<3> nghbr_point (static_cast<float*> (data_n));
     arma::fvec::fixed<3> initial_point (static_cast<float*> (data_p));
     arma::fvec::fixed<3> initial_normal (static_cast<float*> (normals_+3*point));
     float residual = fabs(arma::dot(initial_normal,(initial_point - nghbr_point)));
         if (residual > residual_threshold_)
           is_a_seed = false;
    }
    return (true);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template <typename M> void
RegionGrowingRGB<M>::getSegmentFromPoint (int index, arma::uvec& cluster)
{
    cluster.clear ();
    bool segmentation_is_possible = initCompute ();
    if ( !segmentation_is_possible )
    {
        deinitCompute ();
        return;
    }

    // first of all we need to find out if this point belongs to cloud
    bool point_was_found = false;
    int number_of_points = static_cast <int> (indices_.size ());
    for (int point = 0; point < number_of_points; point++)
    if ( (indices_)(point) == index)
     {
       point_was_found = true;
       break;
     }

   if (point_was_found)
   {
      if (clusters_.empty ())
     {
       clusters_.clear ();
       point_neighbours_.clear ();
       point_labels_.clear ();
       num_pts_in_segment_.clear ();
       point_distances_.clear ();
       segment_neighbours_.clear ();
       segment_distances_.clear ();
       segment_labels_.clear ();
       number_of_segments_ = 0;

       segmentation_is_possible = prepareForSegmentation ();
       if ( !segmentation_is_possible )
       {
         deinitCompute ();
         return;
       }

       findPointNeighbours ();
       applySmoothRegionGrowingAlgorithm ();
       RegionGrowing<M>::assembleRegions ();

       findSegmentNeighbours ();
       applyRegionMergingAlgorithm ();
     }
     // if we have already made the segmentation, then find the segment
     // to which this point belongs
     std::vector<arma::uvec>::iterator i_segment;
     std::vector<arma::uword> cluster_vec;
     for (i_segment = clusters_.begin (); i_segment != clusters_.end (); i_segment++)
     {
       bool segment_was_found = false;
       for (size_t i_point = 0; i_point < i_segment->size (); i_point++)
       {
         if ((*i_segment)[i_point] == index)
         {
           segment_was_found = true;
           cluster_vec.clear ();
           cluster_vec.reserve(i_segment->size ());
           std::copy (i_segment->begin(), i_segment->end (), std::back_inserter (cluster_vec));
           cluster = arma::uvec(cluster_vec);
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

}
#endif
