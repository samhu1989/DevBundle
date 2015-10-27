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

ColorArray::RGB32 COMMONSHARED_EXPORT ColorArray::DefaultColor[DefaultColorNum_] = {
    {0XFF66CCCC},
    {0XFFCCFF66},
    {0XFFFF99CC},
    {0XFFFF9999},
    {0XFFFFCC99},
    {0XFFFF6666},
    {0XFFFFFF66},
    {0XFF99CC66},
    {0XFF666699},
    {0XFFFF9999},
    {0XFFFF6600},
    {0XFF99CC33},
    {0XFFCC3399},
    {0XFFFF9900},
    {0XFFFF9966},
    {0XFFFF0033},
    {0XFFFFCC00},
    {0XFFFF9900},
    {0XFF99CC33},
    {0XFF993366},
    {0XFFFFFF66},
    {0XFF666633},
    {0XFF66CCCC},
    {0XFF666699}
};
