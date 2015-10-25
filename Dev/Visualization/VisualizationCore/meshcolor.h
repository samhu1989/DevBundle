#ifndef MESHCOLOR_H
#define MESHCOLOR_H
#include <iostream>
#include <visualizationcore_global.h>
namespace ColorArray
{
    extern "C"{
        typedef union{
            struct{
                uint8_t r0;
                uint8_t g0;
                uint8_t b0;
                uint8_t a0;
                uint8_t r1;
                uint8_t g1;
                uint8_t b1;
                uint8_t a1;
            }rgba;
            uint64_t color;
        }RGB64;

        typedef union{
            uint32_t color;
            struct{
                uint8_t b;
                uint8_t g;
                uint8_t r;
                uint8_t a;
            }rgba;
        }RGB32;

        const int32_t  DefaultColorNum_ = 24;
        extern RGB32 DefaultColor[DefaultColorNum_];
    }

    typedef struct RGBArray{
        const uint8_t dim_ = 3;
        long size_ = 0;
        uint8_t* data_=NULL;
        void reset(long size,uint8_t r,uint8_t g, uint8_t b);
        ~RGBArray(){if(data_)delete[]data_;}
        }RGBArray;
}

template <typename M>
class MeshColor
{
    typedef M Mesh;
public:
    MeshColor(const Mesh&);
    MeshColor(const MeshColor &);
    ColorArray::RGBArray vertex_colors_array(void);
    void* vertex_colors(void);
protected:
    ColorArray::RGBArray v_colors;
private:
    const Mesh& Ref_;//the ref Mesh that bound with this custom color
};
#include "MeshColor.hpp"
#endif // MESHCOLOR_H
