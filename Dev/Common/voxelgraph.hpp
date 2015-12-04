#include "voxelgraph.h"
#include "KDtree.hpp"
#include "nanoflann.hpp"
#include "MeshColor.h"
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
    if(voxel_centers.n_cols!=voxel_size.size())return false;
    if(voxel_centers.n_cols!=voxel_colors.n_cols)return false;
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

template <typename M>
void VoxelGraph<M>::sv2pix(arma::Col<uint32_t>&sv,arma::Col<uint32_t>&pix)
{
    arma::uvec vl = voxel_label;
    vl -= arma::uvec(vl.size(),arma::fill::ones);
    arma::uvec out_of_bounds = arma::find( vl >= sv.size() || vl < 0 );
    if(!out_of_bounds.is_empty())std::logic_error("!out_of_bounds.is_empty()");
    if(pix.size()!=vl.size())std::logic_error("pix.size()!=vl.size()");
    pix = sv(vl);
}

using namespace nanoflann;
template <typename M>
void VoxelGraph<M>::match(
        M&mesh,
        std::vector<float>&gscore,
        std::vector<float>&cscore,
        arma::vec&score,
        double dist_th,
        double color_var
        )
{
    MeshKDTreeInterface<M> points(mesh);
    KDTreeSingleIndexAdaptor<
            L2_Simple_Adaptor<float,MeshKDTreeInterface<M>>,
            MeshKDTreeInterface<M>,
            3,arma::uword>
            kdtree(3,points,KDTreeSingleIndexAdaptorParams(5));
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
    double min_dist = std::numeric_limits<double>::max();
    for( size_t p_i = 0 ; p_i < Ref_.n_vertices() ; ++ p_i )
    {
        kdtree.knnSearch(&pts[3*p_i],5,search_idx.memptr(),search_dist.memptr());
        arma::fvec current_c;
        ColorArray::RGB2Lab(ref_c_mat.col(p_i),current_c);
        arma::fmat nearest_c;
        ColorArray::RGB2Lab(m_c_mat.cols(search_idx),nearest_c);
        nearest_c.each_col() -= current_c;
        arma::frowvec color_dist = arma::sqrt(arma::sum(arma::square(nearest_c)));
        arma::uword min_idx;
        color_dist.min(min_idx);
        arma::uword sv_idx = voxel_label(p_i) - 1;
        arma::uword match_idx = search_idx(min_idx);
        if(match_idx>gscore.size())std::logic_error("match_idx>gscore.size()");
        if(match_idx>cscore.size())std::logic_error("match_idx>cscore.size()");
        if( search_dist(min_idx) < dist_th )
        {
            sv_match_score(sv_idx) += 1.0/(1.0+search_dist(min_idx))/(1.0+(color_dist(min_idx)/color_var));
        }
        if( min_dist > search_dist(min_idx) ) min_dist = search_dist(min_idx);
        sv_geo_score(sv_idx) += gscore[match_idx];
        sv_color_score(sv_idx) += cscore[match_idx];
    }
    if( min_dist > dist_th )std::cerr<<"min_dist:"<<min_dist<<std::endl;
    for(size_t sv_i = 0 ; sv_i < voxel_centers.n_cols ; ++ sv_i )
    {
        sv_match_score(sv_i) /= voxel_size(sv_i);
        sv_geo_score(sv_i) /= voxel_size(sv_i);
        sv_color_score(sv_i) /= voxel_size(sv_i);
    }
//    double match_max = arma::max(sv_match_score);
//    if(match_max!=0.0)sv_match_score /= match_max;
//    else {
//        std::cerr<<"All zeros in sv_match_score with dist th="<<dist_th<<std::endl;
//    }
//    sv_geo_score /= arma::max(sv_geo_score);
//    sv_color_score /= arma::max(sv_color_score);
    if(!sv_match_score.is_finite())std::cerr<<"Infinite in sv_match_score"<<std::endl;
    for(size_t sv_i = 0 ; sv_i < voxel_centers.n_cols ; ++ sv_i )
    {
        score(sv_i) = sv_match_score(sv_i)*sv_geo_score(sv_i)*sv_color_score(sv_i);
    }
}

template <typename M>
double VoxelGraph<M>::voxel_similarity(size_t v1,size_t v2,double dist_th,double color_var)
{
    if(v1>=voxel_centers.n_cols)std::logic_error("v1>=voxel_centers.n_cols");
    if(v2>=voxel_centers.n_cols)std::logic_error("v2>=voxel_centers.n_cols");
    //use the distance between closest point as distance
    double spatial_dist = std::numeric_limits<double>::max();
    arma::fmat pts((float*)Ref_.points(),3,Ref_.n_vertices(),false,true);
    arma::uvec idx1 = arma::find( voxel_label == (v1+1) );
    arma::uvec idx2 = arma::find( voxel_label == (v2+1) );
    if(idx1.is_empty())std::logic_error("idx1.is_empty()");
    if(idx2.is_empty())std::logic_error("idx2.is_empty()");
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
    double color_dist;
    arma::fvec Lab1;
    ColorArray::RGB2Lab(voxel_colors.col(v1),Lab1);
    arma::fvec Lab2;
    ColorArray::RGB2Lab(voxel_colors.col(v2),Lab2);
    color_dist = arma::norm( Lab1 - Lab2 );
    return std::exp( -( color_dist / color_var ) )*std::exp( -spatial_dist );
}

template <typename M>
void VoxelGraph<M>::get_XYZLab(arma::fmat&voxels,const arma::uvec&indices)
{
    arma::fmat XYZ;
    arma::fmat Lab;
    if(!indices.is_empty())
    {
        XYZ = voxel_centers.cols(indices);
        ColorArray::RGB2Lab(voxel_colors.cols(indices),Lab);
    }else{
        XYZ = voxel_centers;
        ColorArray::RGB2Lab(voxel_colors,Lab);
    }
    voxels = arma::join_cols(XYZ,Lab);
}

template <typename M>
void VoxelGraph<M>::get_Lab(arma::fmat&voxels,const arma::uvec&indices)
{
    arma::fmat Lab;
    if(!indices.is_empty())
    {
        ColorArray::RGB2Lab(voxel_colors.cols(indices),Lab);
    }else{
        ColorArray::RGB2Lab(voxel_colors,Lab);
    }
    voxels = Lab;
}
