#include "voxelgraph.h"
#include "KDtree.hpp"
#include "nanoflann.hpp"

template <typename M>
bool VoxelGraph<M>::save(const std::string&path)
{
    if(!voxel_centers.save(path+"\\centers.fvec.arma",arma::arma_binary))return false;
    if(!voxel_colors.save(path+"\\colors.Mat_uint8_t.arma",arma::arma_binary))return false;
    if(!voxel_size.save(path+"\\sizes.uvec.arma",arma::arma_binary))return false;
    if(!voxel_neighbors.save(path+"\\neighbors.Mat_uint16_t.arma",arma::arma_binary))return false;
    if(!voxel_label.save(path+"\\labels.uvec.arma",arma::arma_binary))return false;
    return true;
}

template <typename M>
bool VoxelGraph<M>::load(const std::string&path)
{
    if(!voxel_centers.load(path+"\\centers.fvec.arma"))return false;
    if(!voxel_colors.load(path+"\\colors.Mat_uint8_t.arma"))return false;
    if(!voxel_size.load(path+"\\sizes.uvec.arma"))return false;
    if(!voxel_neighbors.load(path+"\\neighbors.Mat_uint16_t.arma"))return false;
    if(!voxel_label.load(path+"\\labels.uvec.arma"))return false;
    return true;
}

template <typename M>
void VoxelGraph<M>::sv2pix(arma::uvec& sv,arma::uvec& pix)
{
    if(sv.size()!=voxel_centers.n_cols){
        std::cerr<<"Can't translate a supervoxel label that was not based on this graph"<<std::endl;
        std::logic_error("sv.size()!=voxel_centers.n_cols");
    }
    if(pix.size()!=voxel_label.size())pix = arma::uvec(voxel_label.size(),arma::fill::zeros);
    else pix.fill(0);
    for(int l = 1 ; l <= voxel_centers.n_cols ; ++l )
    {
        arma::uvec indices = arma::find(voxel_label==l);
        pix(indices).fill( sv(l-1) );
    }
}

using namespace nanoflann;
template <typename M>
void VoxelGraph<M>::match(
        M&mesh,
        std::vector<float>&gscore,
        std::vector<float>&cscore,
        arma::vec&score
        )
{
    MeshKDTreeInterface<M> points(mesh);
    KDTreeSingleIndexAdaptor<
            L2_Simple_Adaptor<float,MeshKDTreeInterface<M>>,
            MeshKDTreeInterface<M>,
            3,arma::uword>
            kdtree(3,points,KDTreeSingleIndexAdaptorParams(9));
    kdtree.buildIndex();
    size_t sv_N = voxel_centers.n_cols;
    arma::uvec search_idx(5);
    arma::fvec search_dist(5);

    arma::vec sv_match_score(sv_N,arma::fill::zeros);
    arma::vec sv_geo_score(sv_N,arma::fill::zeros);
    arma::vec sv_color_score(sv_N,arma::fill::zeros);
    score.resize(sv_N);

    float* pts = (float*)Ref_.points();
    arma::Mat<uint8_t> ref_c_mat((uint8_t*)Ref_.vertex_colors(),3,Ref_.n_vertices(),false,true);
    arma::Mat<uint8_t> m_c_mat((uint8_t*)mesh.vertex_colors(),3,mesh.n_vertices(),false,true);
    for( size_t p_i = 0 ; p_i < Ref_.n_vertices() ; ++ p_i )
    {
        kdtree.knnSearch(&pts[3*p_i],5,search_idx.memptr(),search_dist.memptr());
        arma::fvec current_c = arma::conv_to<arma::fvec>::from(ref_c_mat.col(p_i));
        arma::fmat nearest_c = arma::conv_to<arma::fmat>::from(m_c_mat.cols(search_idx));
        nearest_c.each_col() -= current_c;
        arma::frowvec color_dist = arma::sum(arma::square(nearest_c));
        arma::uword min_idx;
        color_dist.min(min_idx);
        arma::uword sv_idx = voxel_label(p_i) - 1;
        if( search_dist(min_idx) < 0.1 )
        sv_match_score(sv_idx) += exp( -search_dist(min_idx) )*exp( -color_dist(min_idx) / 255.0 );
        sv_geo_score(sv_idx) += gscore[p_i];
        sv_color_score(sv_idx) += cscore[p_i];
    }
    for(size_t sv_i = 0 ; sv_i < voxel_centers.n_cols ; ++ sv_i )
    {
        sv_match_score(sv_i) /= voxel_size(sv_i);
        sv_geo_score(sv_i) /= voxel_size(sv_i);
        sv_color_score(sv_i) /= voxel_size(sv_i);
    }
    sv_geo_score /= arma::max(sv_geo_score);
    sv_color_score /= arma::max(sv_color_score);
    for(size_t sv_i = 0 ; sv_i < voxel_centers.n_cols ; ++ sv_i )
    {
        score(sv_i) = sv_match_score(sv_i)*sv_geo_score(sv_i)*sv_color_score(sv_i);
    }
}

template <typename M>
double VoxelGraph<M>::voxel_similarity(size_t v1,size_t v2)
{
    if(v1>=voxel_centers.n_cols)std::logic_error("v1>=voxel_centers.n_cols");
    if(v2>=voxel_centers.n_cols)std::logic_error("v2>=voxel_centers.n_cols");
    //use the distance between closest point as distance
    double spatial_dist = std::numeric_limits<double>::max();
    arma::fmat pts((float*)Ref_.points(),3,Ref_.n_vertices(),false,true);
    arma::uvec idx1 = arma::find( voxel_label == (v1+1) );
    arma::uvec idx2 = arma::find( voxel_label == (v2+1) );
    arma::fmat pts1;
    arma::fmat pts2;
    if( idx1.size() > idx2.size() )
    {
        pts1 = pts.cols(idx1);
        pts2 = pts.cols(idx2);
    }else{
        pts1 = pts.cols(idx2);
        pts2 = pts.cols(idx1);
    }
    ArmaKDTreeInterface<arma::fmat> search_pts(pts1);
    KDTreeSingleIndexAdaptor<
            L2_Simple_Adaptor<float,ArmaKDTreeInterface<arma::fmat>>,
            ArmaKDTreeInterface<arma::fmat>,
            3,arma::uword>
            kdtree(3,search_pts,KDTreeSingleIndexAdaptorParams(1));
    kdtree.buildIndex();
    float* p = (float*)pts2.memptr();
    arma::uword i;
    float d;
    for( size_t pi = 0 ; pi < pts2.n_cols ; ++ pi )
    {
        kdtree.knnSearch(&p[3*pi],1,&i,&d);
        if( d < spatial_dist )spatial_dist = d;
    }
    //use the center color as the color
    double rgb_dist;
    arma::vec rgb1 = arma::conv_to<arma::vec>::from(voxel_colors.col(v1));
    arma::vec rgb2 = arma::conv_to<arma::vec>::from(voxel_colors.col(v2));
    rgb_dist = arma::norm( rgb1 - rgb2 ) / 255.0;
    return exp( - rgb_dist )*exp( - spatial_dist );
}
