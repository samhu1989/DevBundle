#ifndef COMMON_H
#define COMMON_H
#include "common_global.h"
#include "MeshType.h"
#include "MeshColor.h"
#include "KDtree.hpp"
#include "Octree.hpp"
#include "configure.h"
#include "mbb.h"
#include "voxelgraph.h"
#ifndef M_PI
#  define M_PI 3.1415926535897932
#endif
namespace ColorArray {
void COMMONSHARED_EXPORT hsv2rgb(float h,float s,float v,RGB32&rgba);
}
template class COMMONSHARED_EXPORT MeshColor<DefaultMesh>;
template class COMMONSHARED_EXPORT VoxelGraph<DefaultMesh>;
template class COMMONSHARED_EXPORT MeshBundle<DefaultMesh>;
typedef MeshOctreeContainer<DefaultMesh> DefaultOctreeContainer;
template class COMMONSHARED_EXPORT unibn::Octree<arma::fvec,DefaultOctreeContainer>;
typedef unibn::Octree<arma::fvec,DefaultOctreeContainer> DefaultOctree;
void COMMONSHARED_EXPORT getRotationFromTwoUnitVectors(const arma::fvec&,const arma::fvec&,arma::fmat&);
inline void fitPlane(
        arma::fvec &center,
        arma::fmat &neighbor,
        arma::fvec &normal,
        float &curvature,
        float& distToOrigin)
{
    neighbor.each_col() -= center;
    arma::fmat c = arma::cov( neighbor.t() );
    arma::fmat U,V;
    arma::fvec s;
    if(!arma::svd(U,s,V,c,"std"))
    {
        normal.fill(std::numeric_limits<float>::quiet_NaN());
        distToOrigin = std::numeric_limits<float>::quiet_NaN();
        curvature = std::numeric_limits<float>::quiet_NaN();
    }else{
        arma::uword minidx;
        curvature = s.min(minidx);
        normal = U.col(minidx);
        distToOrigin = - arma::dot(normal,center);
    }
}
#endif // COMMON_H
