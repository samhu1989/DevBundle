#include "MeshColor.h"
#include "visualizationcore_global.h"
void VISUALIZATIONCORESHARED_EXPORT ColorArray::RGBArray::reset(long size, uint8_t r, uint8_t g, uint8_t b)
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
