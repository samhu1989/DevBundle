#ifndef COLORHISTOGRAM_HPP
#define COLORHISTOGRAM_HPP
#include "colorhistogram.h"
#include "common.h"
#include <armadillo>
namespace Feature
{
template<typename Mesh>
ColorHistogramLab<Mesh>::ColorHistogramLab(
        uint64_t NL,
        uint64_t Na,
        uint64_t Nb
        ):NL_(NL),Na_(Na),Nb_(Nb)
{
    stepL_ = ( ColorArray::Lab_L_max - ColorArray::Lab_L_min ) / NL_;
    stepa_ = ( ColorArray::Lab_ab_max - ColorArray::Lab_ab_min ) / Na_;
    stepb_ = ( ColorArray::Lab_ab_max - ColorArray::Lab_ab_min ) / Nb_;
}
template<typename Mesh>
void ColorHistogramLab<Mesh>::extract(Mesh&m,arma::vec&hist)
{
    hist = arma::vec( (NL_+ 1)*(Na_+1)*(Nb_+1) ,arma::fill::ones);

    arma::Mat<uint8_t> rgb_mat((uint8_t*)m.vertex_colors(),3,m.n_vertices(),false,true);
    arma::fmat Lab_mat;
    ColorArray::RGB2Lab(rgb_mat,Lab_mat);
    for( size_t index = 0 ; index < Lab_mat.n_cols ; ++index )
    {
        arma::fvec Lab = Lab_mat.col(index);
        hist[( indexL(Lab(0))*Na_ + indexa(Lab(1)) )*Nb_ + indexb(Lab(2))] += 1.0 ;
    }
    hist /= arma::accu(hist);
}
template<typename Mesh>
void ColorHistogramLab<Mesh>::extract(Mesh& m,arma::fvec& feature)
{
    arma::vec ff;
    extract(m,ff);
    feature = arma::conv_to<arma::fvec>::from(ff);
}
template<typename Mesh>
ColorHistogramRGB<Mesh>::ColorHistogramRGB(
        uint64_t Nr,
        uint64_t Ng,
        uint64_t Nb
        ):Nr_(Nr),Ng_(Ng),Nb_(Nb)
{
    stepr_ = 255.0 / Nr_;
    stepg_ = 255.0 / Ng_;
    stepb_ = 255.0 / Nb_;
}
template<typename Mesh>
void ColorHistogramRGB<Mesh>::extract(Mesh&m,arma::vec&hist)
{
    hist = arma::vec( Nr_*Ng_*Nb_ ,arma::fill::ones);
    arma::Mat<uint8_t> rgb_mat((uint8_t*)m.vertex_colors(),3,m.n_vertices(),false,true);
    arma::fmat frgb_mat = arma::conv_to<arma::fmat>::from(rgb_mat);
    for( size_t index = 0 ; index < frgb_mat.n_cols ; ++index )
    {
        arma::fvec rgb = frgb_mat.col(index);
        hist[( indexr(rgb(0))*Ng_ + indexg(rgb(1)) )*Nb_ + indexb(rgb(2))] += 1.0 ;
    }
    hist /= arma::accu(hist);
}
template<typename Mesh>
void ColorHistogramRGB<Mesh>::extract(Mesh& m,arma::fvec& feature)
{
    arma::vec ff;
    extract(m,ff);
    feature = arma::conv_to<arma::fvec>::from(ff);
}
}
#endif

