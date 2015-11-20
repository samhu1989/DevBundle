#ifndef COMMON_H
#define COMMON_H
#include "common_global.h"
#include "MeshType.h"
#include "MeshColor.h"
#include "KDtree.hpp"
#include "Octree.hpp"
#include "configure.h"
#include "mbb.h"
#ifndef M_PI
#  define M_PI 3.1415926535897932
#endif
namespace ColorArray {
void COMMONSHARED_EXPORT hsv2rgb(float h,float s,float v,RGB32&rgba);
}
template class COMMONSHARED_EXPORT MeshColor<DefaultMesh>;
typedef MeshOctreeContainer<DefaultMesh> DefaultOctreeContainer;
template class COMMONSHARED_EXPORT unibn::Octree<arma::fvec,DefaultOctreeContainer>;
typedef unibn::Octree<arma::fvec,DefaultOctreeContainer> DefaultOctree;

#endif // COMMON_H
