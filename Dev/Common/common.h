#ifndef COMMON_H
#define COMMON_H
#include "common_global.h"
#include "MeshType.h"
#include "MeshColor.h"
namespace ColorArray {

void COMMONSHARED_EXPORT hsv2rgb(float h,float s,float v,RGB32&rgba);

}
#endif // COMMON_H
