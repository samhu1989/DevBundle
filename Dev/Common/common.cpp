#include "common.h"
#include "common_global.h"
#include "MeshType.h"
#include "MeshColor.h"


void ColorArray::hsv2rgb(float h,float s,float v,RGB32&rgba)
{
int hi = (int(h) / 60) % 6;
float f = h / 60.0 - hi;
uint8_t p = 255.0*v*(1-s);
uint8_t q = 255.0*v*(1.0-f*s);
uint8_t t = 255.0*v*(1.0-(1.0-f)*s);
uint8_t v_ = 255.0*v;
    switch(hi)
    {
    case 0:
        rgba.rgba.r =v_;rgba.rgba.g=t;rgba.rgba.b = p;return;
    case 1:
        rgba.rgba.r =q;rgba.rgba.g=v_;rgba.rgba.b = p;return;
    case 2:
        rgba.rgba.r=p;rgba.rgba.g=v_;rgba.rgba.b=t;return;
    case 3:
        rgba.rgba.r=p;rgba.rgba.g=q;rgba.rgba.b=v_;return;
    case 4:
        rgba.rgba.r=t;rgba.rgba.g=p;rgba.rgba.b=v_;return;
    case 5:
        rgba.rgba.r=v_;rgba.rgba.g=p;rgba.rgba.b=q;return;
    }
}
