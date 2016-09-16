#ifndef ROBUSTREGION_HPP
#define ROBUSTREGION_HPP
#include "robustregion.h"

namespace Segmentation{
template<typename Mesh>
RobustRegionDetection<Mesh>::RobustRegionDetection():NormalizedCuts<Mesh>()
{
    n_base_segments_ = 5;
}
template<typename Mesh>
void RobustRegionDetection<Mesh>::cutGraph(typename MeshBundle<Mesh>::Ptr m,arma::uvec& label)
{
    computeW_Graph(m);
    generate_base_segments();
    solve_consensus_segment(m,label);
}
template<typename Mesh>
bool RobustRegionDetection<Mesh>::configure(Config::Ptr config)
{
    if(!NormalizedCuts<Mesh>::configure(config))return false;
    if(config->has("Base_Seg_Num")){
        n_base_segments_ = config->getInt("Base_Seg_Num");
    }
    if(!peac.configure(config))return false;
    return true;
}
template<typename Mesh>
void RobustRegionDetection<Mesh>::generate_base_segments(typename MeshBundle<Mesh>::Ptr m)
{
    computeW_Graph(m);
    generate_base_segments();
}

template<typename Mesh>
void RobustRegionDetection<Mesh>::generate_base_segment(const QImage& img)
{
    computeW_Image(img);
    generate_base_segments();
}

template<typename Mesh>
void RobustRegionDetection<Mesh>::solve_consensus_segment(typename MeshBundle<Mesh>::Ptr m, arma::uvec& label)
{
    arma::umat edge = arma::conv_to<arma::umat>::from(m->graph_.voxel_neighbors);
    assert(base_segments_.n_rows==m->graph_.voxel_centers.n_cols);
    peac.compute(base_segments_.t(),edge,label);
}

template<typename Mesh>
void RobustRegionDetection<Mesh>::solve_consensus_segment(const QImage& img,arma::uvec& label)
{
    std::vector<arma::uword> edge_vec;
    arma::uword N = img.height()*img.width();
    edge_vec.reserve(10*N);
    for(int r=0;r<img.height();r++)
    {
//        std::cerr<<"r:"<<r<<std::endl;
        for(int c=0;c<img.width();c++)
        {
            uint32_t wi = r*img.width()+c;
            for(int i=std::max(0,r-5);i<std::min(img.height(),r+5);i++)
            {
                for(int j=std::max(0,c-5);j<std::min(img.width(),c+5);j++)
                {
                    uint32_t wj = i*img.width()+j;
                    if(wi<wj)
                    {
                        edge_vec.push_back(wi);
                        edge_vec.push_back(wj);
                    }
                }
            }
        }
    }
    arma::umat edge(edge_vec.data(),2,edge_vec.size()/2,true,true);
    assert(base_segments_.n_rows==( img.width()*img.height() ));
    peac.compute(base_segments_.t(),edge,label);
}

template<typename Mesh>
void RobustRegionDetection<Mesh>::generate_base_segments()
{
    assert(W_&&(0<W_.use_count()));
    if( base_segments_.n_rows!=W_->n_rows || base_segments_.n_cols!=n_base_segments_)
    {
        base_segments_ = arma::umat(W_->n_rows,n_base_segments_,arma::fill::zeros);
    }
    for(arma::uword i_segments_ = 0 ; i_segments_ < n_base_segments_ ; ++ i_segments_ )
    {
        decomposeGPS();
        clustering_Kmean();
        base_segments_.col(i_segments_) = (gmm_.assign(Y_.t(),arma::eucl_dist)).t();
    }
}
}
#endif // ROBUSTREGION_H
