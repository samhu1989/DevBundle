#include "MeshColor.h"
#include "common_global.h"

void COMMONSHARED_EXPORT ColorArray::RGBArray::reset(long size, uint8_t r, uint8_t g, uint8_t b)
{
    int N;
    if(size_!=size)
    {
        if(data_)delete[] data_;
        size_ = size;
        N = 1 + ( size_/ 2 );
        data_ = (uint8_t*)(new uint64_t[ N ]); // aligned allocation
    }
    RGB64 rgb64;
    rgb64.rgba.r0 = r;
    rgb64.rgba.g0 = g;
    rgb64.rgba.b0 = b;
    rgb64.rgba.a0 = 255;
    rgb64.rgba.r1 = r;
    rgb64.rgba.g1 = g;
    rgb64.rgba.b1 = b;
    rgb64.rgba.a1 = 255;

    uint64_t* long_ptr = (uint64_t*)data_;
    for(long i = 0 ; i < N ; i ++ )
    {
        *long_ptr = rgb64.color;
        long_ptr++;
    }
}

void ColorArray::Lab2RGB(const arma::fvec& Lab, arma::Mat<uint8_t>& rgb)
{

}

ColorArray::RGB32 COMMONSHARED_EXPORT ColorArray::DefaultColor[DefaultColorNum_] = {
    {0XFF000000},
    {0XFFBF80FF},
    {0XFFFF8A80},
    {0XFF808AFF},
    {0XFFFF9F80},
    {0XFFFFAA80},
    {0XFFFFB580},
    {0XFFFFBF80},
    {0XFFFFCA80},
    {0XFFFFD480},
    {0XFFFFDF80},
    {0XFFFFEA80},
    {0XFFFFF480},
    {0XFFFFFF80},
    {0XFFF4FF80},
    {0XFFEAFF80},
    {0XFFDFFF80},
    {0XFFD5FF80},
    {0XFFCAFF80},
    {0XFFBFFF80},
    {0XFFB5FF80},
    {0XFF809FFF},
    {0XFF9FFF80},
    {0XFF95FF80},
    {0XFF8AFF80},
    {0XFF80FF80},
    {0XFF80FF95},
    {0XFF80FF9F},
    {0XFF80FFAA},
    {0XFF80FFB5},
    {0XFF80FFBF},
    {0XFF80FFCA},
    {0XFF80FFD4},
    {0XFF80FF95},
    {0XFF80FFDF},
    {0XFF80FF95},
    {0XFF80FFEA},
    {0XFF80FFF4},
    {0XFF80FF95},
    {0XFF80F4FF},
    {0XFF80EAFF},
    {0XFF80DFFF},
    {0XFF80D4FF},
    {0XFF80CAFF},
    {0XFF80BFFF},
    {0XFF80B5FF},
    {0XFF80EAFF},
    {0XFF80AAFF},
    {0XFF80EAFF},
    {0XFF8095FF},
    {0XFF80FFFF},
    {0XFF80FF8A},
    {0XFF8080FF},
    {0XFF9580FF},
    {0XFF9F80FF},
    {0XFFA280FF},
    {0XFFAC80FF},
    {0XFFFF9580},
    {0XFFAAFF80}
};
