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
void RobustRegionDetection<Mesh>::generate_base_segments(
        typename MeshBundle<Mesh>::Ptr m,
        const arma::uvec& mask,
        bool re_use)
{
    masked_.clear();
    if(!mask.empty())
    {
        arma::uvec indices = arma::find(mask==0);
        if(!indices.empty()&&indices.size()!=mask.size())
        {
            m->graph_.getSvIndex(indices,masked_);
        }
        masked_ -= 1;
    }
    computeW_Graph(m);
    generate_base_segments(re_use);
}


template<typename Mesh>
void RobustRegionDetection<Mesh>::generate_base_segments(typename MeshBundle<Mesh>::Ptr m,bool re_use)
{
    computeW_Graph(m);
    generate_base_segments(re_use);
}

template<typename Mesh>
void RobustRegionDetection<Mesh>::generate_base_segment(const QImage& img,bool re_use)
{
    computeW_Image(img);
    generate_base_segments(re_use);
}

template<typename Mesh>
void RobustRegionDetection<Mesh>::solve_consensus_segment(typename MeshBundle<Mesh>::Ptr m, arma::uvec& label)
{
    //using neighbor and neighbor's neighbor
    std::vector<arma::uword> edge_vec;
    std::vector<std::vector<arma::uword>> nbs(m->graph_.voxel_centers.n_cols);
    arma::Mat<uint16_t>::iterator iter;
    //neighbor
    for(iter=m->graph_.voxel_neighbors.begin();iter!=m->graph_.voxel_neighbors.end();)
    {
        arma::uword i = *iter;
        ++iter;
        arma::uword j = *iter;
        ++iter;
        nbs[i].push_back(j);
        nbs[j].push_back(i);
        edge_vec.push_back(i);
        edge_vec.push_back(j);
    }
    //neighbor's neighbor
    arma::sp_mat nbnb(m->graph_.voxel_centers.n_cols,m->graph_.voxel_centers.n_cols);
    std::vector<arma::uword>::iterator nbiter;
    for(iter=m->graph_.voxel_neighbors.begin();iter!=m->graph_.voxel_neighbors.end();)
    {
        arma::uword i = *iter;
        ++iter;
        arma::uword j = *iter;
        ++iter;
        for(nbiter=nbs[i].begin();nbiter!=nbs[i].end();++nbiter)
        {
            if(*nbiter!=j){
                nbnb(*nbiter,j)=1.0;
                nbnb(j,*nbiter)=1.0;
            }
        }
        for(nbiter=nbs[j].begin();nbiter!=nbs[j].end();++nbiter)
        {
            if(*nbiter!=i){
                nbnb(*nbiter,i)=1.0;
                nbnb(i,*nbiter)=1.0;
            }
        }
    }
    for(arma::sp_mat::iterator iter=nbnb.begin();iter!=nbnb.end();++iter)
    {
        if(*iter==1.0&&iter.row()<iter.col()){
            edge_vec.push_back(iter.row());
            edge_vec.push_back(iter.col());
        }
    }
    arma::umat edge(edge_vec.data(),2,edge_vec.size()/2,true,true);
    assert(base_segments_.n_rows==m->graph_.size());
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
void RobustRegionDetection<Mesh>::generate_base_segments(bool re_use_gps)
{
    assert(W_&&(0<W_.use_count()));
    if( base_segments_.n_rows!=W_->n_rows || base_segments_.n_cols!=n_base_segments_)
    {
        base_segments_ = arma::umat(W_->n_rows,n_base_segments_,arma::fill::zeros);
    }
    if(re_use_gps)decomposeGPS();
    for(arma::uword i_segments_ = 0 ; i_segments_ < n_base_segments_ ; ++ i_segments_ )
    {
        if(!re_use_gps)decomposeGPS();
        clustering_Kmean();
        base_segments_.col(i_segments_) = (gmm_.assign(Y_.t(),arma::eucl_dist)).t();
        if(!masked_.empty())
        {
            arma::uvec l = base_segments_.col(i_segments_);
            l += 1;
            l(masked_).fill(0);
            base_segments_.col(i_segments_) = l;
        }
    }
}
}
#endif // ROBUSTREGION_H
