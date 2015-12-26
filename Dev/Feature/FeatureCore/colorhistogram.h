#ifndef COLORHISTOGRAM_H
#define COLORHISTOGRAM_H
#include "common.h"
namespace Feature
{
template<typename Mesh>
class ColorHistogramLab
{
public:
    ColorHistogramLab(
            uint64_t NL,
            uint64_t Na,
            uint64_t Nb
            );
    void extract(const Mesh&,arma::vec&);
    void extract(const Mesh&,arma::fvec&);
protected:
    inline uint64_t indexL(float L)
    {
        return std::round( (L - ColorArray::Lab_L_min) / stepL_ );
    }
    inline uint64_t indexa(float a)
    {
        return std::round( (a - ColorArray::Lab_ab_min) / stepa_ );
    }
    inline uint64_t indexb(float b)
    {
        return std::round( (b - ColorArray::Lab_ab_min) / stepb_ );
    }
private:
    uint64_t NL_;
    uint64_t Na_;
    uint64_t Nb_;
    float stepL_;
    float stepa_;
    float stepb_;
};
template<typename Mesh>
class ColorHistogramRGB
{
public:
    ColorHistogramRGB(
            uint64_t Nr,
            uint64_t Ng,
            uint64_t Nb
            );
    void extract(const Mesh&,arma::vec&);
    void extract(const Mesh&,arma::fvec&);
protected:
    inline uint64_t indexr(float r)
    {
        uint64_t i = std::floor( ( r - std::numeric_limits<float>::epsilon() ) / stepr_ );
        return i>0?i:0;
    }
    inline uint64_t indexg(float g)
    {
        uint64_t i = std::floor( ( g - std::numeric_limits<float>::epsilon() ) / stepg_ );
        return i>0?i:0;
    }
    inline uint64_t indexb(float b)
    {
        uint64_t i = std::floor( ( b - std::numeric_limits<float>::epsilon() )  / stepb_ );
        return i>0?i:0;
    }
private:
    uint64_t Nr_;
    uint64_t Ng_;
    uint64_t Nb_;
    float stepr_;
    float stepg_;
    float stepb_;
};
}
#include "colorhistogram.hpp"
#endif // COLORHISTOGRAM_H
